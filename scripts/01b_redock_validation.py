#!/usr/bin/env python3
"""
Stage 1b: redocking validation.

Removes the co-crystallised ligand, docks it back into its own site, and asks
whether the top-ranked pose reproduces the crystallographic one.

    python scripts/01b_redock_validation.py --target nep

This is the first of two validation tests, and the weaker of them. It answers
"can the protocol reproduce a known geometry?" A pass is necessary but nowhere
near sufficient: a protocol can reproduce a crystal pose perfectly and still be
unable to tell an active from an inactive. That is what 01c tests.

The convention is that a top-pose RMSD below 2.0 A counts as success. This is a
community convention, not a law of nature, and it is generous: 2.0 A is roughly
the width of a phenyl ring.

A note on RMSD. Comparing two poses atom by atom is wrong whenever a molecule
has symmetry. The two oxygens of a carboxylate are chemically identical, so a
pose with them swapped is the same pose, but naive atom-index RMSD scores it as
badly wrong. This script uses RDKit's GetBestRMS, which searches substructure
matches and returns the best. For carboxylate-bearing NEP inhibitors, the
difference between the two methods is routinely more than 1 A, which is the
difference between passing and failing.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    DATA_INTERIM, RESULTS, banner, get_target, load_config, require_binary,
    setup_logging, write_provenance,
)

RMSD_PASS_THRESHOLD = 2.0


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------

def ligand_pdb_to_sdf(ligand_pdb: Path, out_sdf: Path, logger) -> Path:
    """
    Convert the extracted crystal ligand to SDF with hydrogens and bond orders.

    PDB files carry no bond orders, so Open Babel infers them from geometry.
    This is usually right and occasionally not, which is why the script prints
    the perceived SMILES for you to eyeball. If it looks wrong, supply a
    correct SDF yourself with --reference-sdf.
    """
    require_binary("obabel", "conda install -c conda-forge openbabel")
    cmd = ["obabel", str(ligand_pdb), "-O", str(out_sdf), "-h", "--partialcharge",
           "gasteiger"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_sdf.exists():
        raise RuntimeError(f"Open Babel failed:\n{res.stderr[:800]}")
    return out_sdf


def sdf_to_pdbqt(sdf: Path, out_pdbqt: Path) -> Path:
    """Prepare a docking-ready ligand PDBQT, preferring Meeko over Open Babel."""
    import shutil as _sh
    if _sh.which("mk_prepare_ligand.py"):
        cmd = ["mk_prepare_ligand.py", "-i", str(sdf), "-o", str(out_pdbqt)]
    else:
        cmd = ["obabel", str(sdf), "-O", str(out_pdbqt), "-p", "7.4"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_pdbqt.exists():
        raise RuntimeError(f"Ligand PDBQT preparation failed:\n{res.stderr[:800]}")
    return out_pdbqt


# ---------------------------------------------------------------------------
# Docking
# ---------------------------------------------------------------------------

def run_vina(receptor: Path, ligand: Path, grid: dict, out: Path,
             exhaustiveness: int, num_modes: int, seed: int, logger) -> Path:
    require_binary("vina", "conda install -c conda-forge vina")
    cx, cy, cz = grid["center"]
    sx, sy, sz = grid["size"]
    cmd = [
        "vina",
        "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--out", str(out),
    ]
    logger.info(f"Docking with Vina (exhaustiveness={exhaustiveness}, seed={seed})")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Vina failed:\n{res.stderr[:1200]}")
    (out.parent / f"{out.stem}_vina.log").write_text(res.stdout)
    return out


def run_autodock4_note(logger) -> None:
    logger.warning(
        "This target is configured for AutoDock4 with AD4Zn. Zinc pseudo-atom "
        "setup and grid map generation are a multi-step process that depends on "
        "your local AutoDock4 installation; see docs/stage1.md for the exact "
        "sequence. Falling back to Vina for this validation run, which will "
        "under-represent zinc coordination. Treat the RMSD as indicative only."
    )


# ---------------------------------------------------------------------------
# RMSD
# ---------------------------------------------------------------------------

def parse_vina_poses(pdbqt_out: Path, logger) -> list[tuple[float, str]]:
    """
    Split a Vina output PDBQT into (affinity, pdb_block) pairs, ranked.

    Vina writes all modes to one file, separated by MODEL/ENDMDL, with the
    affinity in a REMARK line.
    """
    poses: list[tuple[float, str]] = []
    affinity = None
    block: list[str] = []
    for line in pdbqt_out.read_text().splitlines():
        if line.startswith("MODEL"):
            block, affinity = [], None
        elif line.startswith("REMARK VINA RESULT"):
            try:
                affinity = float(line.split()[3])
            except (IndexError, ValueError):
                affinity = None
        elif line.startswith("ENDMDL"):
            if block:
                poses.append((affinity if affinity is not None else 0.0,
                              "\n".join(block)))
        elif line.startswith(("ATOM", "HETATM")):
            # Strip PDBQT-specific trailing columns to get valid PDB
            poses_line = line[:66]
            block.append(poses_line)
    logger.info(f"Parsed {len(poses)} poses from {pdbqt_out.name}")
    return poses


def best_rms(ref_mol, probe_mol) -> float:
    """
    Symmetry-corrected RMSD.

    GetBestRMS enumerates substructure matches between the two molecules and
    returns the lowest RMSD over all of them, which is what handles equivalent
    atoms (carboxylate oxygens, phenyl ring flips, amide nitrogens).
    """
    from rdkit.Chem import AllChem, rdMolAlign
    ref = AllChem.RemoveHs(ref_mol)
    probe = AllChem.RemoveHs(probe_mol)
    return rdMolAlign.GetBestRMS(probe, ref)


def naive_rms(ref_mol, probe_mol) -> float:
    """Atom-index RMSD, computed only to show how much symmetry matters."""
    import numpy as np
    from rdkit.Chem import AllChem
    ref = AllChem.RemoveHs(ref_mol)
    probe = AllChem.RemoveHs(probe_mol)
    if ref.GetNumAtoms() != probe.GetNumAtoms():
        return float("nan")
    a = ref.GetConformer().GetPositions()
    b = probe.GetConformer().GetPositions()
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def validate(cfg: dict, key: str, logger, exhaustiveness: int | None = None) -> dict:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    block = get_target(cfg, key)
    outdir = DATA_INTERIM / key
    grid_json = outdir / "grid.json"
    receptor = outdir / f"{key}_receptor.pdbqt"
    ligand_pdb = outdir / f"{key}_reference_ligand.pdb"

    for p in (grid_json, receptor, ligand_pdb):
        if not p.exists():
            raise RuntimeError(
                f"{p.name} not found. Run 01a_prepare_targets.py --target {key} first."
            )

    banner(logger, f"REDOCKING VALIDATION  |  {key.upper()}")
    grid = json.loads(grid_json.read_text())

    if block["docking"].get("engine") == "autodock4":
        run_autodock4_note(logger)

    ref_sdf = outdir / f"{key}_reference_ligand.sdf"
    ligand_pdb_to_sdf(ligand_pdb, ref_sdf, logger)

    ref_mol = Chem.MolFromMolFile(str(ref_sdf), removeHs=False)
    if ref_mol is None:
        raise RuntimeError(
            f"RDKit could not read {ref_sdf.name}. Bond perception from the PDB "
            f"likely failed. Supply a correct SDF with --reference-sdf."
        )
    smiles = Chem.MolToSmiles(Chem.RemoveHs(ref_mol))
    logger.info(f"Reference ligand perceived as: {smiles}")
    logger.info("Check that SMILES looks right before trusting anything below.")

    lig_pdbqt = outdir / f"{key}_reference_ligand.pdbqt"
    sdf_to_pdbqt(ref_sdf, lig_pdbqt)

    docked = outdir / f"{key}_redock_out.pdbqt"
    ex = exhaustiveness or block["docking"].get("exhaustiveness", 16)
    run_vina(receptor, lig_pdbqt, grid, docked, ex,
             block["docking"].get("num_modes", 9),
             cfg["project"].get("random_seed", 42), logger)

    poses = parse_vina_poses(docked, logger)
    if not poses:
        raise RuntimeError("Vina produced no poses. Check the grid box covers the site.")

    rows = []
    for rank, (affinity, block_text) in enumerate(poses, start=1):
        probe = Chem.MolFromPDBBlock(block_text, removeHs=False, sanitize=False)
        if probe is None:
            continue
        try:
            rms = best_rms(ref_mol, probe)
            naive = naive_rms(ref_mol, probe)
        except Exception as exc:
            logger.debug(f"Pose {rank}: RMSD failed ({exc})")
            continue
        rows.append({"rank": rank, "affinity_kcal_mol": affinity,
                     "rmsd_A": round(rms, 3), "naive_rmsd_A": round(naive, 3)})

    if not rows:
        raise RuntimeError("No pose could be compared to the reference.")

    top = rows[0]
    best = min(rows, key=lambda r: r["rmsd_A"])
    passed = top["rmsd_A"] <= RMSD_PASS_THRESHOLD

    logger.info("")
    logger.info(f"{'rank':>5}  {'affinity':>9}  {'RMSD':>7}  {'naive RMSD':>11}")
    for r in rows:
        logger.info(f"{r['rank']:>5}  {r['affinity_kcal_mol']:>9.2f}  "
                    f"{r['rmsd_A']:>7.2f}  {r['naive_rmsd_A']:>11.2f}")
    logger.info("")
    logger.info(f"Top pose RMSD      : {top['rmsd_A']:.2f} A")
    logger.info(f"Best pose RMSD     : {best['rmsd_A']:.2f} A (rank {best['rank']})")
    logger.info(f"Threshold          : {RMSD_PASS_THRESHOLD:.1f} A")
    logger.info(f"RESULT             : {'PASS' if passed else 'FAIL'}")

    if not passed:
        logger.error(
            "Redocking failed. Do not proceed to screening. Usual causes, in "
            "order of likelihood: grid box does not cover the site or is too "
            "small; wrong ligand HET code; the catalytic metal was dropped "
            "during receptor preparation; bond orders misperceived from the PDB."
        )
    if best["rank"] != 1 and best["rmsd_A"] <= RMSD_PASS_THRESHOLD:
        logger.warning(
            f"The correct pose was found but ranked {best['rank']}, not first. "
            f"The sampling works; the scoring function is the weak link. Expect "
            f"the same failure mode across the screen."
        )

    report = {
        "target": key, "pdb_id": grid.get("pdb_id"),
        "reference_smiles": smiles,
        "exhaustiveness": ex,
        "top_pose_rmsd_A": top["rmsd_A"],
        "best_pose_rmsd_A": best["rmsd_A"],
        "best_pose_rank": best["rank"],
        "threshold_A": RMSD_PASS_THRESHOLD,
        "passed": passed,
        "poses": rows,
    }
    out_json = RESULTS / f"validation_redock_{key}.json"
    out_json.write_text(json.dumps(report, indent=2))
    write_provenance(out_json, "01b_redock_validation",
                     {"target": key, "exhaustiveness": ex, "grid": grid},
                     tools=["vina", "obabel"])
    logger.info(f"Report written: {out_json.relative_to(RESULTS.parent)}")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", default="all")
    ap.add_argument("--exhaustiveness", type=int, default=None,
                    help="Override config; higher is slower and more thorough")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logger = setup_logging("01b", args.verbose)
    cfg = load_config()
    keys = list(cfg["targets"]) if args.target == "all" else [args.target]

    results, failed = [], []
    for key in keys:
        try:
            results.append(validate(cfg, key, logger, args.exhaustiveness))
        except Exception as exc:
            logger.error(f"{key}: {exc}")
            failed.append(key)

    print()
    for r in results:
        logger.info(f"{r['target']:>6}: {'PASS' if r['passed'] else 'FAIL'} "
                    f"(top pose {r['top_pose_rmsd_A']:.2f} A)")
    if failed or any(not r["passed"] for r in results):
        return 1
    logger.info("Next:  python scripts/01c_enrichment_benchmark.py --target nep")
    return 0


if __name__ == "__main__":
    sys.exit(main())
