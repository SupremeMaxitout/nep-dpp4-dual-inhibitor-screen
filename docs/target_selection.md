# Choosing and verifying target structures

`config/targets.yaml` ships with candidate PDB entries marked `# VERIFY`. They are
starting points, not validated choices. This document explains how to confirm or
replace them.

The original project's exact accession codes were not preserved, so this step has to
be redone rather than reproduced.

## Criteria

A structure is suitable for this pipeline if it meets all of the following:

1. **Human protein.** Rodent orthologues differ enough in the S1 pocket to change
   the ranking.
2. **Holo, not apo.** A co-crystallised inhibitor is required. It defines the grid
   box and provides the redocking reference. An apo structure gives you neither, and
   its side chains are often in a conformation no ligand can occupy.
3. **Resolution 2.5 A or better.** Below this, side-chain positions in the site are
   not reliable enough to dock against.
4. **Complete binding site.** Check for missing residues or disordered loops near
   the site. The REMARK 465 records list them.
5. **For neprilysin: the catalytic zinc must be present.** An NEP structure without
   its zinc is useless here.

## Procedure

1. Search RCSB by UniProt accession: P08473 for NEP, P27487 for DPP-IV, P12821 for ACE.
2. Filter to X-ray, human, resolution better than 2.5 A, with a bound ligand.
3. Open the candidate entry and record: the ligand HET code, the chain identifier, the
   resolution, and for NEP the zinc-coordinating residue numbers as that entry numbers
   them. Residue numbering is not consistent across entries.
4. Enter these into `config/targets.yaml`, replacing every `VERIFY` placeholder.
5. Run `python scripts/01b_redock_validation.py --target nep`. If redocking RMSD comes
   back above 2.0 A, the structure or the grid box is wrong. Do not proceed past a
   failing validation.

## Recording the choice

Whatever you settle on, write it down in `docs/target_selection.md` under a short
"chosen structures" heading, with the date and the reason. A reader should be able to
see not only which structure you used but why you rejected the alternatives. This is
the difference between a reproducible pipeline and a set of scripts.
