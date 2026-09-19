# Dual Neprilysin / DPP-IV Inhibitor Discovery Pipeline

A reproducible structure-based virtual screening pipeline for identifying molecules
with dual inhibitory activity against **neprilysin (NEP)** and **dipeptidyl
peptidase-IV (DPP-IV)**, as a candidate strategy for type 2 diabetes in patients at
risk of heart failure with reduced ejection fraction (HFrEF).

---

## Why this target pair

Sacubitril/valsartan established neprilysin inhibition as a mainstay of HFrEF therapy,
while DPP-IV inhibition is a well-tolerated glucose-lowering mechanism. Patients sit at
the intersection of both conditions frequently, and a single molecule acting on both
targets would reduce pill burden and avoid the drug-drug interaction surface of
combination therapy.

The pharmacological rationale is attractive. The medicinal chemistry is not
straightforward, and this repository is as much a documentation of *why* as a
documentation of *how*. See [`docs/methods_rationale.md`](docs/methods_rationale.md).

---

## Provenance and scope

**Please read this section before citing or evaluating this work.**

This pipeline originates from my final-year B.Pharm project (University of Mumbai),
carried out several years before this repository was created. The original scripts and
intermediate data were not preserved under version control.

This repository is therefore a **reimplementation**, not an archive. It reconstructs the
original methodology from my records and recollection, and in several places it
deliberately **improves on what was actually done at the time**. Every such divergence
is documented explicitly in
[`docs/differences_from_original.md`](docs/differences_from_original.md) so that the
historical work and the present implementation are never confused with one another.

**Attribution.** In the original project, the machine-learning model training and the
generative expansion of scaffold libraries were performed by a PhD colleague, not by me.
My contribution covered target preparation, docking validation, database mining and
library construction, the docking and rescoring cascade, the cross-target selection
logic, and scaffold analysis. This repository reflects my portion of the work. The ML
stage is described for completeness but is not reimplemented here.

---

## What the pipeline does

```
  [1] Target preparation and docking validation
        Redock co-crystallised ligands, confirm pose recovery,
        benchmark the protocol against property-matched decoys
                              |
  [2] Library construction
        Mine ChEMBL and ZINC for annotated actives per target,
        filter on potency, molecular weight and lipophilicity
                              |
  [3] Primary docking and rescoring
        Dock each library against its own target, select poses,
        rescore, rank by both affinity and ligand efficiency
                              |
  [4] Cross-docking and dual-activity selection
        Dock each surviving set against the opposite target,
        keep only molecules active against both
                              |
  [5] Scaffold analysis and constrained enumeration
        Extract shared scaffolds, enumerate analogues under an
        explicit dual-pharmacophore constraint, counter-screen
        against ACE and DPP-8/9
                              |
              Prioritised set for in vitro testing
```

Each stage is a self-contained, independently runnable script with its own
documentation. Nothing is a black box.

---

## Repository layout

| Path | Contents |
|------|----------|
| `scripts/` | Numbered, independently runnable pipeline stages |
| `config/` | Target definitions, grid boxes, filter thresholds |
| `docs/` | Method rationale, stage-by-stage explanation, divergences from the original |
| `data/raw/` | Downloaded structures and database dumps (git-ignored) |
| `data/interim/` | Intermediate artefacts: PDBQT files, docking logs |
| `data/processed/` | Final tabulated results |
| `results/` | Figures, ranked tables, reports |
| `notebooks/` | Optional exploratory analysis |

---

## Installation

The pipeline depends on both Python packages and standalone binaries. The conda
environment handles the former; the latter are described in
[`docs/installation.md`](docs/installation.md).

```bash
git clone https://github.com/<YOUR-USERNAME>/nep-dpp4-dual-inhibitor-screen.git
cd nep-dpp4-dual-inhibitor-screen
conda env create -f environment.yml
conda activate nepdpp4
python scripts/00_check_environment.py
```

`00_check_environment.py` verifies every dependency and tells you exactly what is
missing and where to get it. Run it first.

---

## Running the pipeline

Every stage script accepts `--limit N` to truncate the working set.

```bash
# Quick demonstration run, finishes on a laptop
python scripts/03_primary_docking.py --target nep --limit 500

# Full-scale run, assumes a cluster
python scripts/03_primary_docking.py --target nep
```

Full-scale execution across both libraries is on the order of 10^5 docking
calculations per target and is intended for HPC. SLURM submission templates are
provided under `scripts/slurm/`.

---

## Compute requirements

| Stage | Demo (`--limit 500`) | Full scale |
|-------|----------------------|------------|
| 1. Validation | ~15 min, laptop | ~1 h, laptop |
| 2. Library build | ~20 min, network-bound | ~2 h, network-bound |
| 3. Primary docking | ~1 h, 8 cores | ~2,000 CPU-hours |
| 4. Cross-docking | ~30 min, 8 cores | ~500 CPU-hours |
| 5. Enumeration and counter-screen | ~1 h, 8 cores | ~800 CPU-hours |

---

## Limitations

Stated plainly, because a virtual screen that does not state its limitations should not
be trusted:

- No compound from this pipeline has been synthesised or assayed. All results are
  computational predictions.
- Docking scores correlate weakly with experimental binding affinity. They are used
  here for triage and ranking, not for quantitative prediction.
- The dual-pharmacophore requirement (an anionic zinc-binding group for NEP, a basic
  nitrogen for DPP-IV) implies zwitterionic character, which carries a well-known
  permeability penalty. This is a chemical constraint the pipeline surfaces rather
  than solves.
- Protein flexibility is not modelled during docking.

---

## Citation

If you use or refer to this pipeline, please cite it via the metadata in
[`CITATION.cff`](CITATION.cff), or use the "Cite this repository" button on GitHub.

---

## Licence

MIT. See [`LICENSE`](LICENSE).

Third-party tools invoked by this pipeline (AutoDock Vina, AutoDock4, Open Babel,
RDKit, and the ChEMBL and ZINC databases) carry their own licences and citation
requirements. These are listed in [`docs/installation.md`](docs/installation.md) and
should be cited alongside this work.
