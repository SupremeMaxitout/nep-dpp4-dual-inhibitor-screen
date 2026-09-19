# Installation

This pipeline uses a mix of Python packages and standalone scientific binaries. The
conda environment covers the former automatically. The latter need a few manual steps,
explained below with the reasoning for each choice.

Run `python scripts/00_check_environment.py` at any point to see what is present and
what is missing.

---

## 1. Conda environment

If you do not already have conda, install
[Miniforge](https://github.com/conda-forge/miniforge) rather than Anaconda. It is
smaller, community-governed, and defaults to the conda-forge channel that this
environment depends on.

```bash
conda env create -f environment.yml
conda activate nepdpp4
```

This gives you RDKit, Open Babel, AutoDock Vina, Meeko, and the data-handling stack.
For the demonstration-scale runs, this is sufficient on its own.

---

## 2. AutoDock4 and the zinc parameter set

**Needed for:** neprilysin docking only.

### Why this is not optional

Neprilysin is a zinc metallopeptidase. Its catalytic Zn²⁺ ion is coordinated by
His583, His587 and Glu646, and essentially every potent inhibitor binds by donating a
ligand to that zinc: a carboxylate in sacubitrilat, a thiol in thiorphan, a
hydroxamate in others.

AutoDock Vina's scoring function has no dedicated term for metal coordination. It sees
the zinc as a generic atom with van der Waals and electrostatic character. The result
is a systematic mis-ranking of exactly the interaction that determines potency at this
target.

AutoDock4 with the **AD4Zn** force field addresses this by placing four tetrahedrally
arranged pseudo-atoms around the zinc, which reproduce its directional coordination
preference. This is why the pipeline uses two docking engines rather than one: Vina for
DPP-IV, where it performs well, and AutoDock4/AD4Zn for neprilysin, where Vina does not.

### Installation

```bash
# Linux / macOS
wget https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-x86_64Linux2.tar
tar -xf autodocksuite-4.2.6-x86_64Linux2.tar
sudo mv x86_64Linux2/autodock4 x86_64Linux2/autogrid4 /usr/local/bin/
```

The AD4Zn parameter file (`AD4Zn.dat`) ships with the suite. The pipeline expects it at
a path you set in `config/targets.yaml`; `00_check_environment.py` will report where it
found it.

> If you are on Apple Silicon, AutoDock4 is x86 only and needs Rosetta 2, or you can
> run the whole pipeline in the provided container.

---

## 3. ADFR Suite (for receptor preparation)

**Needed for:** converting cleaned PDB files into PDBQT receptors.

Download from <https://ccsb.scripps.edu/adfr/downloads/>. The pipeline uses
`prepare_receptor` from this suite.

**A note on why not AutoDockTools.** The original project used the legacy
`prepare_receptor4.py` and `prepare_ligand4.py` scripts from AutoDockTools, which
depend on Python 2 and MGLTools. Those are effectively unmaintained. This
reimplementation uses `prepare_receptor` from ADFR for proteins and **Meeko** for
ligands, both of which are actively maintained, Python 3 native, and produce equivalent
output. If you have a working MGLTools installation, the legacy path still works; see
`scripts/01_prepare_targets.py --legacy`.

---

## 4. Pose selection and rescoring (IGModel)

**Needed for:** stage 3 pose selection and affinity rescoring.

The original project used IGModel's RMSD branch for pose selection and its pKd branch
for affinity rescoring. Availability and installation for this tool should be verified
against its current upstream repository before you rely on it, as its distribution has
changed over time.

The pipeline therefore treats rescoring as a **pluggable step**. `config/targets.yaml`
exposes a `rescoring.backend` key with these options:

| Backend | Requirement | Notes |
|---------|-------------|-------|
| `vina_native` | none | Vina's own affinity. Baseline; always available. |
| `igmodel` | external install | Matches the original methodology. |
| `rfscore_vs` | ODDT (`pip install oddt`) | Machine-learned rescoring, well validated |
| `mmgbsa` | AmberTools or OpenMM | Slowest, most physical. Top hits only. |

Nothing downstream depends on which you choose. The default is `vina_native` so that a
fresh clone runs end to end without any external setup.

---

## 5. Optional: containerised run

A `Dockerfile` is provided for reproducibility and for machines where the binaries above
are awkward to install.

```bash
docker build -t nepdpp4 .
docker run -it -v $(pwd):/work nepdpp4
```

---

## Citing the tools

If you publish anything derived from this pipeline, these need citing alongside it:

- **AutoDock Vina 1.2** — Eberhard J, Santos-Martins D, Tillack AF, Forli S.
  *J Chem Inf Model* 2021.
- **AutoDock4 / AD4Zn** — Santos-Martins D, Forli S, Ramos MJ, Olson AJ.
  *J Chem Inf Model* 2014.
- **Open Babel** — O'Boyle NM et al. *J Cheminform* 2011.
- **RDKit** — Landrum G, open-source cheminformatics, <https://www.rdkit.org>.
- **Meeko** — Forli Lab, Scripps Research.
- **ChEMBL** — Zdrazil B et al. *Nucleic Acids Res* 2024.
- **ZINC** — Irwin JJ et al. *J Chem Inf Model* 2020.

Verify these against current records before submission; version and year details change.
