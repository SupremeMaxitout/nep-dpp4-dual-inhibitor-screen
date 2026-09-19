# Stage 1: target preparation and docking validation

Three scripts, run in order. Nothing downstream is meaningful unless all three pass.

```bash
python scripts/01a_prepare_targets.py --target all
python scripts/01b_redock_validation.py --target nep
python scripts/01c_enrichment_benchmark.py --target nep --limit 100
```

---

## Before you start: resolve the config placeholders

`01a` refuses to run while any `VERIFY` marker remains in `config/targets.yaml`.
This is deliberate and fatal rather than a warning, because a pipeline that runs
happily against the wrong structure produces output that looks completely normal
and is worthless.

See [`target_selection.md`](target_selection.md) for how to confirm the PDB
entries, chains, ligand HET codes and zinc-coordinating residue numbers.

---

## 01a: target preparation

Downloads the PDB entry, keeps one chain, discards waters and crystallisation
additives, **retains functional metal ions**, splits out the co-crystallised
ligand, derives the grid box from that ligand's centroid, and converts the
receptor to PDBQT.

### The zinc

Neprilysin's catalytic Zn²⁺ is the single most important atom in the site. A
generic "remove all HETATM records" cleanup deletes it, and every downstream
number is then wrong in a way that is invisible.

The script keeps a whitelist of functional ions (`FUNCTIONAL_HETATMS`) and, for
any target with `metal.present: true`, checks after PDBQT conversion that the
metal actually survived. If it did not, you get a loud error rather than a
silently broken receptor.

### The grid box

Centre comes from the reference ligand's centroid. Size comes from the config,
but is grown if the ligand's own bounding box plus 4 Å padding is larger. A box
smaller than the reference ligand cannot reproduce the crystal pose, which would
make the next step meaningless.

Derived values are written to `data/interim/<target>/grid.json`. If you want to
override the centre, set it explicitly in the config and the script will use
yours and record that it did.

### Outputs

```
data/interim/nep/
├── nep_receptor.pdb          cleaned protein + zinc
├── nep_receptor.pdbqt        docking-ready receptor
├── nep_reference_ligand.pdb  crystal pose, for validation
├── grid.json                 box definition + provenance
└── grid.json.meta.json       config, git commit, tool versions
```

---

## 01b: redocking validation

Takes the crystal ligand out, docks it back, asks whether the top pose
reproduces the original geometry. Pass threshold is 2.0 Å RMSD.

### On RMSD and symmetry

Comparing two poses atom-by-atom is wrong whenever the molecule has symmetry.
The two oxygens of a carboxylate are chemically indistinguishable, so a pose
with them swapped **is the same pose**, but naive index-based RMSD scores it as
badly wrong.

Tested on a phenylacetic acid where only the two carboxylate oxygens are
exchanged:

| Method | RMSD |
|--------|------|
| Naive atom-index | 0.995 Å |
| Symmetry-corrected (`GetBestRMS`) | 0.000 Å |

For the carboxylate-bearing chemotypes this project depends on, that difference
is routinely the difference between a reported pass and a reported fail. The
script uses the symmetry-corrected value and reports the naive one alongside, so
you can see the effect on your own ligand.

### Reading the result

**Top pose under 2.0 Å.** Pass. Proceed.

**Correct pose found but ranked second or lower.** The script warns explicitly.
Sampling works; the scoring function is the weak link. Expect the same failure
across the whole screen, and plan to rescore.

**All poses above 2.0 Å.** Stop. In rough order of likelihood: the grid box does
not cover the site or is too small; the ligand HET code is wrong; the catalytic
metal was dropped during preparation; Open Babel misperceived bond orders from
the PDB. The script prints the perceived SMILES for exactly this reason. Look at
it before you look at anything else.

### AutoDock4 and zinc

`config/targets.yaml` sets neprilysin's engine to `autodock4` for AD4Zn support.
The zinc pseudo-atom workflow is installation-specific, so `01b` currently warns
and falls back to Vina for the validation run, marking the result as indicative.
Full AD4Zn integration lands with Stage 3, where it matters for ranking rather
than for geometry.

---

## 01c: enrichment benchmark

The test the original project omitted, and the one that decides whether the
screen is worth running.

### Why redocking is not enough

Redocking shows the protocol can reproduce a geometry it was effectively handed
the answer to. It says nothing about whether the protocol can rank an active
above an inactive, which is the only thing a virtual screen actually does.

### Property-matched decoys

Mix known actives with compounds chosen to look physicochemically identical but
be presumed inactive, dock everything, and see whether the actives rise.

Each decoy must match an active within tolerance on molecular weight (±25 Da),
logP (±1.0), HBD (±1), HBA (±2), rotatable bonds (±2) and formal charge (exact),
while sharing no more than 0.35 Tanimoto similarity with any active.

The matching constraint is what makes the test honest. Random decoys would
differ from the actives in weight and charge, and docking would separate them on
those trivial grounds. You would measure beautiful enrichment and learn nothing.

The dissimilarity constraint stops a genuine active sneaking in as a decoy,
which would poison the benchmark in the other direction.

### Interpreting the numbers

| ROC-AUC | Verdict | What to do |
|---------|---------|------------|
| < 0.60 | Unusable | Screening with this is an expensive random draw. Fix the protocol. |
| 0.60–0.70 | Weak | Coarse triage only. Don't trust ranking within the survivors. |
| 0.70–0.85 | Acceptable | A normal docking result. Proceed. |
| > 0.85 | Verify | Unusually good. Check for a property leak in the decoys first. |

**EF1%** matters more than AUC in practice, because you only ever buy or
synthesise from the top of the list. EF1% of 10 means the top 1% holds ten times
more actives than chance would give. Its ceiling is the total/actives ratio, so
with 50 actives among 1000 compounds the best possible EF1% is 20.

**BEDROC (α=20)** weights the top ~8% of the list. Reported because AUC rewards
good behaviour in the tail, which nobody cares about.

### Threshold calibration

Rather than inheriting the original project's fixed −9.0 kcal/mol cutoff, the
script reports the score capturing the configured percentile of your own ranked
list, and writes it to `results/validation_enrichment_<target>.json` for Stage 3
to pick up. A cutoff derived from your protocol's own behaviour beats a number
copied from someone else's paper.

### Getting a background set

You need a broad, property-diverse set of presumed inactives. A random subset of
ZINC's in-stock catalogue works. Save it as `data/raw/background.csv` with a
`smiles` column. Around 50,000 compounds gives comfortable matching; fewer and
the script will warn you that the decoy set is thin.

Actives come from Stage 2, or from any CSV you supply with a `smiles` column.

---

## Provenance

Every artefact gets a sibling `.meta.json` recording the config that produced it,
the git commit, timestamps and tool versions. This costs nothing and is the only
reason anyone, including you, will be able to defend these results in a year.

This repository exists because that record was not kept the first time.
