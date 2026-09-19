#!/usr/bin/env python3
"""
Stage 1a: target preparation.

Turns a raw PDB entry into something dockable, and derives the search box from
the co-crystallised ligand rather than from guesswork.

    python scripts/01a_prepare_targets.py --target nep
    python scripts/01a_prepare_targets.py --target all

What it does, in order:

  1. Download the PDB entry (cached in data/raw/).
  2. Keep one chain. Discard waters and crystallisation additives.
  3. KEEP the catalytic zinc if the target has one. This is the step that is
     most often got wrong: a generic "strip all heteroatoms" cleanup removes
     the single most important feature of the neprilysin site.
  4. Split out the co-crystallised ligand into its own file. It serves two
     purposes downstream: the redocking reference pose, and the grid centre.
  5. Compute the grid box from the ligand's centroid.
  6. Convert the receptor to PDBQT.

Outputs land in data/interim/<target>/ and the grid definition is written to
data/interim/<target>/grid.json for the docking stages to read.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    DATA_RAW, DATA_INTERIM, banner, check_no_placeholders, get_target,
    load_config, require_binary, setup_logging, write_provenance,
)

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Ions and cofactors that are part of the protein's function and must survive
# cleanup. Everything else in HETATM records is assumed to be crystallisation
# debris (buffers, cryoprotectants, detergents) and is discarded.
FUNCTIONAL_HETATMS = {"ZN", "MG", "CA", "MN", "FE", "NA", "CU", "CO", "NI"}


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def fetch_pdb(pdb_id: str, logger) -> Path:
    pdb_id = pdb_id.upper()
    dest = DATA_RAW / f"{pdb_id}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        logger.info(f"Using cached structure {dest.name}")
        return dest
    url = RCSB_URL.format(pdb_id=pdb_id)
    logger.info(f"Downloading {pdb_id} from RCSB")
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            dest.write_bytes(resp.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"RCSB returned HTTP {exc.code} for '{pdb_id}'. Check the accession "
            f"code in config/targets.yaml against https://www.rcsb.org/"
        ) from exc
    logger.info(f"Saved {dest.name} ({dest.stat().st_size / 1024:.0f} kB)")
    return dest


# ---------------------------------------------------------------------------
# Selection classes for Biopython's PDBIO
# ---------------------------------------------------------------------------

class ReceptorSelect(Select):
    """Protein of one chain, plus functional metal ions. No waters, no ligand."""

    def __init__(self, chain_id: str, keep_hetatms: set[str]):
        self.chain_id = chain_id
        self.keep_hetatms = keep_hetatms

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        hetflag, _, _ = residue.id
        if hetflag == " ":
            return True                                    # standard amino acid
        if hetflag == "W" or residue.get_resname().strip() == "HOH":
            return False                                   # water
        return residue.get_resname().strip() in self.keep_hetatms

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


class LigandSelect(Select):
    """The co-crystallised reference ligand only."""

    def __init__(self, chain_id: str, hetcode: str):
        self.chain_id = chain_id
        self.hetcode = hetcode.strip().upper()

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.get_resname().strip().upper() == self.hetcode

    def accept_atom(self, atom):
        return atom.get_altloc() in (" ", "A")


# ---------------------------------------------------------------------------
# Grid box
# ---------------------------------------------------------------------------

def grid_from_ligand(ligand_pdb: Path, size, padding: float = 4.0) -> dict:
    """
    Derive a docking box from the reference ligand.

    Centre is the ligand centroid. Size is whatever the config specifies, but
    if the ligand's own bounding box plus padding is larger, the box is grown
    to fit. A box smaller than the reference ligand cannot reproduce the
    crystal pose, which makes validation meaningless.
    """
    coords = []
    with open(ligand_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]),
                               float(line[46:54])])
    if not coords:
        raise RuntimeError(f"No atoms found in {ligand_pdb}")

    arr = np.asarray(coords)
    centre = arr.mean(axis=0)
    span = arr.max(axis=0) - arr.min(axis=0) + 2 * padding
    final = [max(float(s), float(sp)) for s, sp in zip(size, span)]

    return {
        "center": [round(float(c), 3) for c in centre],
        "size": [round(f, 1) for f in final],
        "reference_ligand_atoms": len(coords),
        "ligand_bounding_box": [round(float(s), 2) for s in span],
    }


# ---------------------------------------------------------------------------
# PDBQT conversion
# ---------------------------------------------------------------------------

def receptor_to_pdbqt(receptor_pdb: Path, out_pdbqt: Path, logger) -> Path:
    """
    Convert a cleaned receptor to PDBQT.

    Preferred route is ADFR's prepare_receptor. Open Babel is a fallback that
    works but handles protonation less carefully, so it warns.
    """
    import shutil as _sh

    if _sh.which("prepare_receptor"):
        cmd = ["prepare_receptor", "-r", str(receptor_pdb),
               "-o", str(out_pdbqt), "-A", "hydrogens", "-U", "nphs_lps"]
        logger.info("Preparing receptor with ADFR prepare_receptor")
    else:
        require_binary("obabel", "conda install -c conda-forge openbabel")
        logger.warning(
            "prepare_receptor not found, falling back to Open Babel. This works, "
            "but ADFR handles histidine protonation and charge assignment more "
            "carefully. See docs/installation.md section 3."
        )
        cmd = ["obabel", str(receptor_pdb), "-O", str(out_pdbqt),
               "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(
            f"Receptor preparation failed.\nCommand: {' '.join(cmd)}\n"
            f"{result.stderr[:1000]}"
        )
    return out_pdbqt


def verify_metal_retained(pdbqt: Path, element: str, logger) -> bool:
    """
    Confirm the catalytic metal survived conversion.

    Some preparation routes silently drop metal ions. For neprilysin that
    invalidates everything downstream, so it is checked rather than assumed.
    """
    text = pdbqt.read_text()
    found = any(
        line.startswith(("ATOM", "HETATM")) and element.upper() in line[76:78].upper()
        for line in text.splitlines()
    )
    if found:
        logger.info(f"Catalytic {element} retained in receptor PDBQT")
    else:
        logger.error(
            f"Catalytic {element} is ABSENT from {pdbqt.name}. Docking against "
            f"this receptor will mis-rank every zinc-binding inhibitor. Check "
            f"the HET code and chain in config/targets.yaml."
        )
    return found


# ---------------------------------------------------------------------------
# Main per-target routine
# ---------------------------------------------------------------------------

def prepare_one(cfg: dict, key: str, logger, force: bool = False) -> dict:
    block = get_target(cfg, key)
    check_no_placeholders(block, key)

    struct_cfg = block["structure"]
    pdb_id = struct_cfg["pdb_id"].upper()
    chain_id = struct_cfg["chain"]
    hetcode = struct_cfg["reference_ligand_hetcode"].strip().upper()

    banner(logger, f"{key.upper()}  |  {block['full_name']}  |  PDB {pdb_id}")

    outdir = DATA_INTERIM / key
    outdir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = outdir / f"{key}_receptor.pdb"
    ligand_pdb = outdir / f"{key}_reference_ligand.pdb"
    receptor_pdbqt = outdir / f"{key}_receptor.pdbqt"
    grid_json = outdir / "grid.json"

    if receptor_pdbqt.exists() and not force:
        logger.info("Already prepared. Use --force to redo.")
        return json.loads(grid_json.read_text())

    raw = fetch_pdb(pdb_id, logger)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(raw))
    model = next(structure.get_models())          # first model only

    chains = [c.id for c in model]
    if chain_id not in chains:
        raise RuntimeError(
            f"Chain '{chain_id}' not present in {pdb_id}. Available: {chains}"
        )

    # Confirm the reference ligand is actually there before doing any work
    hetcodes = {
        r.get_resname().strip().upper()
        for r in model[chain_id]
        if r.id[0] not in (" ", "W")
    }
    if hetcode not in hetcodes:
        raise RuntimeError(
            f"Reference ligand '{hetcode}' not found in chain {chain_id} of "
            f"{pdb_id}.\nHeteroatoms present: {sorted(hetcodes)}\n"
            f"Fix reference_ligand_hetcode in config/targets.yaml."
        )

    keep = set(FUNCTIONAL_HETATMS)
    metal_cfg = block.get("metal", {}) or {}
    has_metal = bool(metal_cfg.get("present"))
    metal_element = (metal_cfg.get("element") or "ZN").upper()

    io = PDBIO()
    io.set_structure(structure)

    io.save(str(receptor_pdb), ReceptorSelect(chain_id, keep))
    logger.info(f"Receptor written: {receptor_pdb.name}")

    io.save(str(ligand_pdb), LigandSelect(chain_id, hetcode))
    logger.info(f"Reference ligand written: {ligand_pdb.name}")

    grid = grid_from_ligand(
        ligand_pdb,
        block["docking"]["grid"]["size"],
        padding=4.0,
    )
    configured_centre = block["docking"]["grid"].get("center")
    if configured_centre:
        logger.info(f"Overriding derived centre with config value {configured_centre}")
        grid["center"] = configured_centre
        grid["center_source"] = "config"
    else:
        grid["center_source"] = "derived from reference ligand centroid"

    grid["spacing"] = block["docking"]["grid"].get("spacing", 0.375)
    grid["pdb_id"] = pdb_id
    grid["chain"] = chain_id
    grid["reference_ligand_hetcode"] = hetcode

    grid_json.write_text(json.dumps(grid, indent=2))
    logger.info(
        f"Grid box: centre {grid['center']}  size {grid['size']} A "
        f"({grid['center_source']})"
    )

    receptor_to_pdbqt(receptor_pdb, receptor_pdbqt, logger)
    logger.info(f"Receptor PDBQT written: {receptor_pdbqt.name}")

    if has_metal:
        grid["metal_retained"] = verify_metal_retained(
            receptor_pdbqt, metal_element, logger
        )
        grid_json.write_text(json.dumps(grid, indent=2))
        if block["docking"].get("engine") == "autodock4":
            logger.info(
                "This target docks with AutoDock4 + AD4Zn. The zinc pseudo-atom "
                "step runs in 01b before docking; see docs/stage1.md."
            )

    write_provenance(
        grid_json, stage="01a_prepare_targets",
        params={"target": key, "pdb_id": pdb_id, "chain": chain_id,
                "hetcode": hetcode, "grid": grid},
        tools=["obabel", "prepare_receptor"],
    )
    return grid


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="all",
                    help="Target key from config/targets.yaml, or 'all'")
    ap.add_argument("--force", action="store_true",
                    help="Re-prepare even if outputs already exist")
    ap.add_argument("--include-antitargets", action="store_true",
                    help="With --target all, also prepare ACE and DPP-8")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logger = setup_logging("01a", args.verbose)
    cfg = load_config()

    if args.target == "all":
        keys = list(cfg["targets"])
        if args.include_antitargets:
            keys += list(cfg.get("antitargets", {}))
    else:
        keys = [args.target]

    failures = []
    for key in keys:
        try:
            prepare_one(cfg, key, logger, force=args.force)
        except Exception as exc:
            logger.error(f"{key}: {exc}")
            failures.append(key)

    print()
    if failures:
        logger.error(f"Failed: {', '.join(failures)}")
        return 1
    logger.info("All targets prepared.")
    logger.info("Next:  python scripts/01b_redock_validation.py --target nep")
    return 0


if __name__ == "__main__":
    sys.exit(main())
