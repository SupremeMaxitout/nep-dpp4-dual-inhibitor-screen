#!/usr/bin/env python3
"""
Stage 0: environment check.

Run this first. It verifies every dependency the pipeline needs, reports what is
missing, and tells you where to get it. Nothing else in the pipeline will give you
a useful error message if a binary is absent, so start here.

    python scripts/00_check_environment.py

Exit code 0 means the pipeline can run at demonstration scale. Optional components
are reported but do not cause failure.
"""

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# What we expect to find.
#
# Each entry: (import name, human label, why it is needed, required?)
# --------------------------------------------------------------------------

PYTHON_PACKAGES = [
    ("rdkit", "RDKit", "molecule handling, fingerprints, scaffolds, SMARTS", True),
    ("pandas", "pandas", "tabular data throughout the pipeline", True),
    ("numpy", "NumPy", "numerical work", True),
    ("yaml", "PyYAML", "reading config/targets.yaml", True),
    ("tqdm", "tqdm", "progress bars on long runs", True),
    ("sklearn", "scikit-learn", "ROC-AUC and enrichment metrics in stage 1", True),
    ("matplotlib", "matplotlib", "figures", True),
    ("Bio", "Biopython", "PDB parsing and cleaning", True),
    ("chembl_webresource_client", "ChEMBL client", "database mining in stage 2", True),
    ("meeko", "Meeko", "ligand PDBQT preparation", True),
    ("MDAnalysis", "MDAnalysis", "trajectory and structure utilities", False),
    ("oddt", "ODDT", "optional RF-Score-VS rescoring backend", False),
]

# Each entry: (executable, human label, why, required?, where to get it)
BINARIES = [
    ("vina", "AutoDock Vina", "primary docking engine (DPP-IV)", True,
     "conda install -c conda-forge vina"),
    ("obabel", "Open Babel", "chemical file format conversion", True,
     "conda install -c conda-forge openbabel"),
    ("autodock4", "AutoDock4", "zinc-aware docking for neprilysin", False,
     "see docs/installation.md section 2"),
    ("autogrid4", "AutoGrid4", "grid map generation for AutoDock4", False,
     "see docs/installation.md section 2"),
    ("prepare_receptor", "ADFR prepare_receptor", "receptor PDBQT preparation", False,
     "see docs/installation.md section 3"),
]

GREEN, RED, YELLOW, BOLD, RESET = (
    "\033[92m", "\033[91m", "\033[93m", "\033[1m", "\033[0m"
)


def tick(ok: bool, required: bool) -> str:
    if ok:
        return f"{GREEN}  found{RESET}"
    return f"{RED}MISSING{RESET}" if required else f"{YELLOW}absent{RESET}"


def check_python_packages() -> list:
    print(f"\n{BOLD}Python packages{RESET}")
    print("-" * 78)
    missing = []
    for module, label, why, required in PYTHON_PACKAGES:
        try:
            importlib.import_module(module)
            ok = True
        except ImportError:
            ok = False
            if required:
                missing.append(label)
        flag = "" if required else " (optional)"
        print(f"  [{tick(ok, required)}]  {label:<22} {why}{flag}")
    return missing


def check_binaries() -> list:
    print(f"\n{BOLD}Standalone binaries{RESET}")
    print("-" * 78)
    missing = []
    for exe, label, why, required, howto in BINARIES:
        ok = shutil.which(exe) is not None
        if not ok and required:
            missing.append((label, howto))
        flag = "" if required else " (optional)"
        print(f"  [{tick(ok, required)}]  {label:<22} {why}{flag}")
        if not ok:
            print(f"                 -> {howto}")
    return missing


def check_layout() -> None:
    print(f"\n{BOLD}Repository layout{RESET}")
    print("-" * 78)
    root = Path(__file__).resolve().parent.parent
    for rel in ["config/targets.yaml", "data/raw", "data/interim",
                "data/processed", "results", "docs"]:
        path = root / rel
        print(f"  [{tick(path.exists(), True)}]  {rel}")
        if not path.exists() and not path.suffix:
            path.mkdir(parents=True, exist_ok=True)
            print(f"                 -> created")


def report_vina_version() -> None:
    if shutil.which("vina") is None:
        return
    try:
        out = subprocess.run(["vina", "--version"], capture_output=True,
                             text=True, timeout=10).stdout.strip()
        print(f"\n  Vina reports: {out.splitlines()[0] if out else 'unknown'}")
    except Exception:
        pass


def main() -> int:
    print(f"{BOLD}NEP / DPP-IV dual-inhibitor pipeline: environment check{RESET}")
    print(f"Python {sys.version.split()[0]} at {sys.executable}")

    missing_pkgs = check_python_packages()
    missing_bins = check_binaries()
    check_layout()
    report_vina_version()

    print("\n" + "=" * 78)
    if not missing_pkgs and not missing_bins:
        print(f"{GREEN}{BOLD}Environment is ready.{RESET}")
        print("\nNext step:  python scripts/01_prepare_targets.py --help")
        return 0

    print(f"{RED}{BOLD}Required components are missing.{RESET}\n")
    if missing_pkgs:
        print("  Python packages: " + ", ".join(missing_pkgs))
        print("  Fix:  conda env create -f environment.yml && conda activate nepdpp4\n")
    for label, howto in missing_bins:
        print(f"  {label}\n  Fix:  {howto}\n")
    print("Items marked 'absent' rather than 'MISSING' are optional. The pipeline")
    print("runs without them at reduced fidelity; see docs/installation.md.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
