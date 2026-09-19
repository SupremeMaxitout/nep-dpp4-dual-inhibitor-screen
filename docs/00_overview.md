# Pipeline overview

This document explains what each stage does, why it exists, and what it hands to the
next stage. Read it before running anything. Each stage also has its own detailed
document, linked below.

The design principle throughout: **every stage produces a human-readable table, and
every filter can be justified.** A virtual screen that cannot explain why a compound
survived is not a screen, it is a lottery.

---

## Stage 0 — Environment check

`scripts/00_check_environment.py`

Verifies dependencies and reports what is missing. Not part of the science; part of
making the science runnable by someone else.

---

## Stage 1 — Target preparation and docking validation

`scripts/01_prepare_targets.py`, `scripts/01b_validate_docking.py`

**What it does.** Downloads the target structures, strips waters and crystallisation
additives while retaining the catalytic zinc in neprilysin, adds hydrogens at
physiological pH, converts to PDBQT, and derives the grid box from the co-crystallised
ligand.

Then it validates in two ways:

1. **Redocking.** The co-crystallised ligand is removed and docked back. If the top
   pose reproduces the crystallographic pose to within 2.0 A RMSD, the protocol can
   at least reproduce a known answer.

2. **Enrichment.** A set of known actives is docked alongside property-matched decoys
   (compounds resembling the actives physicochemically but presumed inactive). The
   protocol is scored by ROC-AUC and enrichment factor at 1%.

**Why the second step matters.** Redocking tests geometry. Enrichment tests
discrimination. A protocol can reproduce a crystal pose perfectly and still rank
actives no better than chance. The original project did the first and not the second,
which is the single most common omission in student virtual screening work.

**Output.** Prepared receptors, grid definitions, a validation report with ROC curves,
and a calibrated affinity threshold to use downstream.

---

## Stage 2 — Library construction

`scripts/02_build_library.py`

**What it does.** Queries ChEMBL for compounds with measured potency against each
target, filters on activity value, assay confidence and physicochemical properties,
strips salts, removes duplicates by InChIKey, discards PAINS, and optionally augments
from ZINC. Produces one library per target.

**Why the property filters are loose.** A strict rule-of-five gate would discard
exactly the zwitterionic chemotypes this project needs, since a molecule carrying both
a carboxylate and a basic amine is polar by construction. The filters are set to
exclude the clearly undevelopable without prejudging the chemistry.

**Output.** Two curated SDF libraries with an accompanying properties table.

---

## Stage 3 — Primary docking and rescoring

`scripts/03_primary_docking.py`, `scripts/03b_rescore.py`

**What it does.** Docks each library against its own target. Neprilysin uses
AutoDock4 with zinc pseudo-atoms; DPP-IV uses Vina. Poses are selected, rescored
through the configured backend, and ranked by both raw affinity and ligand efficiency.

**Why two ranking columns.** Raw docking score scales with molecular size, so a
ranking by score alone is partly a ranking by molecular weight. Ligand efficiency
(score divided by heavy atom count) corrects for this and identifies compounds that
bind efficiently rather than merely largely. Both are reported so you can see the
difference the choice makes.

**Output.** A ranked affinity table per target, and the surviving hit set at the
calibrated threshold.

---

## Stage 4 — Cross-docking and dual-activity selection

`scripts/04_cross_dock.py`

**What it does.** Takes the NEP hit set and docks it against DPP-IV, and vice versa.
Compounds passing the threshold on both targets form the dual-active set.

**Why this replaces a similarity step.** The original project bridged the two
libraries using 2D Tanimoto similarity, then cross-docked only the similar pairs. That
was a reasonable compute-saving heuristic, but fingerprint similarity does not imply a
shared binding mode. Two molecules can be 2D-similar and bind completely differently.
Cross-docking everything answers the actual question directly.

**Output.** The dual-active set, with a per-compound affinity profile across both
targets.

---

## Stage 5 — Scaffold analysis, enumeration, counter-screening

`scripts/05_scaffolds.py`, `scripts/05b_enumerate.py`, `scripts/05c_counterscreen.py`

**What it does.** Extracts Murcko scaffolds from the dual actives and ranks them by
frequency and by binding-site occupancy. Scaffolds are filtered against an explicit
dual-pharmacophore SMARTS requirement: an anionic zinc-binding group for NEP, and a
basic nitrogen for DPP-IV. Surviving scaffolds are enumerated into analogue sets,
filtered on synthetic accessibility, docked, and finally counter-screened against ACE
and DPP-8/9.

**Why the pharmacophore filter is explicit.** In the original project, the hope was
that a molecule satisfying both targets would emerge from the similarity analysis.
Making the requirement an explicit constraint turns a hope into a specification, and
makes the subsequent failure modes interpretable.

**Why the counter-screen is not optional.** Omapatrilat, a dual NEP/ACE inhibitor,
failed clinically because of angioedema. DPP-8/9 cross-reactivity is a known toxicity
liability. Selectivity is part of the answer, not an afterthought.

**Output.** A prioritised, selectivity-filtered compound set with full provenance for
every surviving structure.

---

## What is deliberately absent

**The machine-learning stage.** The original project trained a model to predict
docking score from structure. That work was done by a colleague and is not
reimplemented here. The methodological critique is recorded in
`docs/differences_from_original.md`.

**Molecular dynamics.** MD and MM/GBSA refinement of the final hits would be the
natural next step and are not included. Docking gives a static, rigid-receptor
approximation; anything advanced to synthesis should be checked for pose stability
first.

**Experimental validation.** None. Everything here is a prediction.
