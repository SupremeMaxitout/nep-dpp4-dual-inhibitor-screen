# Differences from the original project

This repository implements a single, improved pipeline. It does not preserve the
original undergraduate implementation verbatim, because the original scripts were not
kept under version control.

This document records every point at which the present implementation departs from what
was actually done, so that the two are never confused. The table is the honest ledger
of the work: what was done then, what is done now, and why the change was made.

---

## Summary table

| # | Stage | Original (B.Pharm project) | This implementation | Reason for change |
|---|-------|---------------------------|---------------------|-------------------|
| 1 | Validation | Redocking of co-crystal ligands, visual pose inspection in UCSF Chimera | Redocking retained, plus retrospective enrichment against property-matched decoys (ROC-AUC, EF1%) | Low redocking RMSD demonstrates geometric reproduction, not the ability to tell actives from inactives. Without an enrichment benchmark there is no evidence the protocol discriminates at all. |
| 2 | NEP docking engine | AutoDock Vina | AutoDock4 with AD4Zn zinc pseudo-atoms | Vina has no metal-coordination term. Neprilysin inhibition is driven by zinc chelation, so the original ranking was systematically unreliable for this target. |
| 3 | File preparation | AutoDockTools / MGLTools (Python 2) | ADFR `prepare_receptor` + Meeko (Python 3) | The legacy toolchain is unmaintained and hard to install. Output is equivalent. |
| 4 | Selection threshold | Fixed cutoff at −9 kcal/mol | Threshold calibrated from the enrichment benchmark, plus a parallel ranking by ligand efficiency | A fixed kcal/mol cutoff is biased toward large, lipophilic molecules, which score well simply for having more atoms. Ligand efficiency normalises for size and yields more developable starting points. |
| 5 | Bridging the two libraries | 2D Tanimoto similarity between NEP hits and DPP-IV hits; similar pairs then cross-docked | Direct cross-docking of both full hit sets against the opposite target; dual actives taken as the intersection | Fingerprint similarity does not imply a shared binding mode. Two molecules can be 2D-similar and bind entirely differently, or be dissimilar and share a pharmacophore. The similarity step was a compute-saving proxy; cross-docking answers the question directly and is now affordable. |
| 6 | Scaffold selection | Scaffolds chosen manually by the project supervisor | Murcko scaffold extraction with frequency and occupancy ranking, plus an explicit dual-pharmacophore SMARTS filter | Makes the selection reproducible rather than expert-dependent, and enforces the chemical requirement that was previously only hoped for (see below). |
| 7 | Analogue enumeration | DataWarrior combinatorial expansion | RDKit-based enumeration under scaffold constraints, with synthetic accessibility scoring | Removes a GUI step from an otherwise scriptable pipeline, and filters out the synthetically implausible structures that unconstrained enumeration produces in bulk. |
| 8 | Antitarget screening | Not performed | ACE and DPP-8/9 counter-screens added | Omapatrilat, a dual NEP/ACE inhibitor, failed clinically due to angioedema. DPP-8/9 cross-reactivity is a known toxicity liability for DPP-IV programmes. For this target pair these counter-screens are not optional. |
| 9 | ML stage | Model trained on docking scores to predict activity from structure (performed by a PhD colleague) | Not reimplemented here | Out of scope for this repository; see note below. |

---

## On the machine-learning stage

In the original project, a model was trained to relate molecular structure to docking
score, in order to prioritise a minimal set of compounds for in vitro testing. That work
was carried out by a PhD colleague and is not my own contribution, so it is not
reimplemented in this repository.

It is worth recording the methodological criticism regardless, since it applies to a
great many published virtual screens:

A model trained on docking scores learns to approximate the **scoring function**, not
the **biology**. Its accuracy ceiling is the accuracy of the docking programme it was
trained against, and docking scores correlate only weakly with experimental affinity. A
model that perfectly reproduced Vina would inherit every one of Vina's errors.

The better approach was available at the time and was overlooked: stage 2 of this
pipeline mines ChEMBL for compounds *filtered by their experimental IC50 values*. Those
experimental potency labels were downloaded and then discarded in favour of computed
scores. A multi-task model trained on experimental pChEMBL values for both targets,
using the docking scores only as an auxiliary signal, would have used data that was
already in hand.

---

## Unrecoverable details

The following details of the original work could not be reconstructed and have been
re-derived rather than reproduced:

- Exact PDB accession codes for the two target structures
- Exact grid box centre coordinates and dimensions
- Exact ChEMBL and ZINC query parameters and retrieval dates
- The identity of the scaffolds selected by the supervisor
- Final ranked compound lists

Where a value had to be chosen afresh, the choice is documented in
`config/targets.yaml` with its justification. Results produced by this repository are
therefore **not** expected to reproduce the original numerical outputs, and should not
be represented as doing so.
