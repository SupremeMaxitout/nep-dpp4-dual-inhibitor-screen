#!/usr/bin/env python3
"""
Stage 1c: enrichment benchmark.

The test the original project never ran, and the one that decides whether the
screen is worth running at all.

    python scripts/01c_enrichment_benchmark.py --target nep --limit 200

The idea
--------
Redocking (01b) shows the protocol can reproduce a geometry it was effectively
given the answer to. It says nothing about whether the protocol can rank an
active above an inactive, which is the only thing a virtual screen actually
does.

So: take compounds known to be active against the target, mix them with decoys
chosen to look physicochemically identical but be presumed inactive, dock the
whole mixture, and see whether the actives float to the top.

Why property-matched decoys
---------------------------
If decoys were drawn at random, they would differ from the actives in molecular
weight, charge and lipophilicity, and docking would separate the two sets on
those trivial grounds. You would measure a beautiful enrichment and learn
nothing. Matching each decoy to an active on molecular weight, logP, charge,
rotatable bonds and hydrogen-bond counts removes that shortcut, while keeping
the topology different so the decoys are unlikely to be real binders.

This is the DUD-E construction principle, applied here directly because DUD-E
has no neprilysin set.

What the numbers mean
---------------------
ROC-AUC     0.5 is random. Below about 0.65 the protocol is not usable.
            0.7 to 0.8 is a typical, workable docking result.
            Above 0.9 on a well-built decoy set is suspicious; check for a
            property leak before celebrating.

EF1%        Enrichment factor in the top 1%. EF1% of 10 means the top 1% of
            the ranked list holds ten times more actives than chance. This
            matters more than AUC for screening, because you only ever buy or
            synthesise from the top of the list.

BEDROC      AUC weighted toward the top of the list, with alpha=20 covering
            roughly the top 8%. Reported because AUC rewards good behaviour
            in the tail, which nobody cares about.

The threshold this stage calibrates
-----------------------------------
Rather than inheriting the original project's fixed -9.0 kcal/mol cutoff, the
script reports the score that captures a chosen percentile of the ranked list,
and writes it back for later stages. A cutoff derived from your own protocol's
behaviour beats a number copied from someone else's paper.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (  # noqa: E402
    DATA_INTERIM, RESULTS, banner, get_target, load_config, setup_logging,
    write_provenance,
)

# Properties on which decoys are matched to actives.
MATCH_PROPERTIES = ["mw", "logp", "hbd", "hba", "rotb", "net_charge"]

# Tolerance for a decoy to count as matching an active on each property.
MATCH_TOLERANCE = {
    "mw": 25.0, "logp": 1.0, "hbd": 1, "hba": 2, "rotb": 2, "net_charge": 0,
}

# Maximum 2D similarity a decoy may share with any active. Above this it may
# genuinely be active, which would poison the benchmark.
MAX_DECOY_SIMILARITY = 0.35


# ---------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------

def compute_properties(smiles: str) -> dict | None:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "smiles": Chem.MolToSmiles(mol),
        "mw": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "net_charge": Chem.GetFormalCharge(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
    }


def fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return gen.GetFingerprint(mol)


# ---------------------------------------------------------------------------
# Decoy selection
# ---------------------------------------------------------------------------

def select_decoys(actives: pd.DataFrame, background: pd.DataFrame,
                  per_active: int, logger) -> pd.DataFrame:
    """
    For each active, pick `per_active` background compounds that match it on
    physicochemical properties but are topologically dissimilar.
    """
    from rdkit import DataStructs

    active_fps = [fingerprint(s) for s in actives["smiles"]]
    active_fps = [f for f in active_fps if f is not None]

    chosen: list[dict] = []
    used: set[str] = set()

    for _, act in actives.iterrows():
        mask = pd.Series(True, index=background.index)
        for prop in MATCH_PROPERTIES:
            tol = MATCH_TOLERANCE[prop]
            mask &= (background[prop] - act[prop]).abs() <= tol
        candidates = background[mask & ~background["smiles"].isin(used)]

        picked = 0
        for _, cand in candidates.iterrows():
            if picked >= per_active:
                break
            fp = fingerprint(cand["smiles"])
            if fp is None:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(fp, active_fps)
            if max(sims) > MAX_DECOY_SIMILARITY:
                continue                       # too similar; could be a real binder
            rec = cand.to_dict()
            rec["matched_active"] = act["smiles"]
            rec["max_similarity_to_actives"] = round(float(max(sims)), 3)
            chosen.append(rec)
            used.add(cand["smiles"])
            picked += 1

        if picked < per_active:
            logger.debug(f"Only {picked}/{per_active} decoys matched one active")

    df = pd.DataFrame(chosen)
    logger.info(f"Selected {len(df)} decoys for {len(actives)} actives "
                f"(target ratio 1:{per_active})")
    if len(df) < len(actives) * per_active * 0.5:
        logger.warning(
            "Fewer than half the requested decoys could be matched. The "
            "background set is probably too small or too narrow in property "
            "space. Enrichment numbers from a thin decoy set are unreliable."
        )
    return df


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def enrichment_factor(labels: np.ndarray, scores: np.ndarray,
                      fraction: float) -> float:
    """
    EF at a given fraction of the ranked list.

    Scores are docking affinities, so more negative is better. The list is
    sorted ascending and the top `fraction` examined.
    """
    n = len(labels)
    n_top = max(1, int(round(n * fraction)))
    order = np.argsort(scores)
    hits_top = labels[order][:n_top].sum()
    hits_all = labels.sum()
    if hits_all == 0:
        return float("nan")
    return float((hits_top / n_top) / (hits_all / n))


def bedroc(labels: np.ndarray, scores: np.ndarray, alpha: float = 20.0) -> float:
    """BEDROC, emphasising early recognition. Truchon & Bayly, 2007."""
    n = len(labels)
    n_act = int(labels.sum())
    if n_act == 0 or n_act == n:
        return float("nan")
    order = np.argsort(scores)
    ranks = np.where(labels[order] == 1)[0] + 1
    ra = n_act / n
    rie_num = np.sum(np.exp(-alpha * ranks / n)) / n_act
    rie_den = (1 / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = rie_num / rie_den
    factor = ra * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
    return float(rie * factor + 1 / (1 - np.exp(alpha * (1 - ra))))


def evaluate(labels: np.ndarray, scores: np.ndarray) -> dict:
    from sklearn.metrics import roc_auc_score
    return {
        "n_total": int(len(labels)),
        "n_actives": int(labels.sum()),
        "n_decoys": int(len(labels) - labels.sum()),
        # Negate: more negative affinity means better, roc_auc expects higher-is-better
        "roc_auc": round(float(roc_auc_score(labels, -scores)), 4),
        "bedroc_alpha20": round(bedroc(labels, scores), 4),
        "ef_0.5pct": round(enrichment_factor(labels, scores, 0.005), 2),
        "ef_1pct": round(enrichment_factor(labels, scores, 0.01), 2),
        "ef_5pct": round(enrichment_factor(labels, scores, 0.05), 2),
    }


def interpret(metrics: dict, logger) -> str:
    auc = metrics["roc_auc"]
    ef1 = metrics["ef_1pct"]
    if auc < 0.6:
        verdict = "UNUSABLE"
        logger.error(
            "The protocol ranks actives barely better than chance. Screening "
            "100k compounds with it would be an expensive random draw. Revisit "
            "the receptor preparation, the grid box, and for neprilysin whether "
            "the zinc is being handled at all."
        )
    elif auc < 0.7:
        verdict = "WEAK"
        logger.warning(
            "Marginal. Usable for coarse triage only. Do not trust the ranking "
            "within the surviving set, and plan to rescore the top hits with "
            "something more physical."
        )
    elif auc < 0.85:
        verdict = "ACCEPTABLE"
        logger.info("A normal, workable docking result. Proceed.")
    else:
        verdict = "STRONG (verify)"
        logger.warning(
            "Unusually good. Before accepting this, check the decoys are truly "
            "property-matched. A leak in molecular weight or charge produces "
            "exactly this number and means nothing."
        )
    logger.info(f"Early enrichment EF1% = {ef1}")
    return verdict


def calibrate_threshold(scores: np.ndarray, percentile: float) -> float:
    """Score at which the top `percentile` of the ranked list is captured."""
    return float(np.percentile(scores, percentile))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True)
    ap.add_argument("--actives", type=Path, default=None,
                    help="CSV with a 'smiles' column. Default: data/raw/<target>_actives.csv")
    ap.add_argument("--background", type=Path, default=None,
                    help="CSV of presumed-inactive compounds to draw decoys from")
    ap.add_argument("--decoys-per-active", type=int, default=50)
    ap.add_argument("--limit", type=int, default=None,
                    help="Cap the number of actives. Use for a quick check.")
    ap.add_argument("--scores", type=Path, default=None,
                    help="Pre-computed docking scores CSV (smiles,score). "
                         "Skips docking and evaluates metrics only.")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logger = setup_logging("01c", args.verbose)
    cfg = load_config()
    key = args.target
    get_target(cfg, key)          # validates the key exists

    banner(logger, f"ENRICHMENT BENCHMARK  |  {key.upper()}")

    actives_path = args.actives or (DATA_INTERIM.parent / "raw" / f"{key}_actives.csv")
    background_path = args.background or (DATA_INTERIM.parent / "raw" / "background.csv")

    if not actives_path.exists():
        logger.error(
            f"No actives file at {actives_path}.\n"
            f"Build one with:  python scripts/02_build_library.py --target {key} "
            f"--potent-only --out {actives_path}\n"
            f"Or supply your own CSV with a 'smiles' column."
        )
        return 1
    if not background_path.exists() and args.scores is None:
        logger.error(
            f"No background set at {background_path}. This should be a broad, "
            f"property-diverse compound set (a random ZINC subset works). See "
            f"docs/stage1.md for how to obtain one."
        )
        return 1

    actives = pd.read_csv(actives_path)
    if "smiles" not in actives.columns:
        logger.error(f"{actives_path.name} has no 'smiles' column")
        return 1
    if args.limit:
        actives = actives.head(args.limit)
        logger.info(f"Limited to {len(actives)} actives")

    props = [compute_properties(s) for s in actives["smiles"]]
    actives = pd.DataFrame([p for p in props if p])
    logger.info(f"{len(actives)} actives with valid structures")

    background = pd.read_csv(background_path)
    bg_props = [compute_properties(s) for s in background["smiles"]]
    background = pd.DataFrame([p for p in bg_props if p])
    logger.info(f"{len(background)} background compounds")

    decoys = select_decoys(actives, background, args.decoys_per_active, logger)

    benchmark = pd.concat([
        actives.assign(label=1, role="active"),
        decoys.assign(label=0, role="decoy"),
    ], ignore_index=True)

    out_csv = DATA_INTERIM / key / f"{key}_benchmark_set.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    benchmark.to_csv(out_csv, index=False)
    logger.info(f"Benchmark set written: {out_csv.name} ({len(benchmark)} compounds)")

    if args.scores is None:
        logger.info("")
        logger.info("Benchmark set built. Dock it, then re-run with --scores:")
        logger.info(f"  python scripts/03_primary_docking.py --input {out_csv} "
                    f"--target {key}")
        logger.info(f"  python scripts/01c_enrichment_benchmark.py --target {key} "
                    f"--scores <docking_results.csv>")
        return 0

    scores_df = pd.read_csv(args.scores)
    merged = benchmark.merge(scores_df[["smiles", "score"]], on="smiles", how="inner")
    logger.info(f"{len(merged)} of {len(benchmark)} compounds have docking scores")
    if len(merged) < len(benchmark) * 0.8:
        logger.warning("More than 20% of the benchmark failed to dock. Metrics "
                       "computed on a biased subset are not trustworthy.")

    labels = merged["label"].to_numpy()
    scores = merged["score"].to_numpy()
    metrics = evaluate(labels, scores)

    logger.info("")
    for k, v in metrics.items():
        logger.info(f"  {k:<16} {v}")
    logger.info("")
    verdict = interpret(metrics, logger)

    pct = cfg["docking"]["selection"].get("target_enrichment_percentile", 1.0)
    threshold = calibrate_threshold(scores, pct)
    logger.info(f"Calibrated threshold at top {pct}%: {threshold:.2f} kcal/mol")
    logger.info(f"(config fallback was "
                f"{cfg['docking']['selection']['affinity_cutoff_kcal_mol']})")

    report = {"target": key, "verdict": verdict, "metrics": metrics,
              "calibrated_threshold_kcal_mol": round(threshold, 2),
              "percentile": pct}
    out_json = RESULTS / f"validation_enrichment_{key}.json"
    out_json.write_text(json.dumps(report, indent=2))
    write_provenance(out_json, "01c_enrichment_benchmark",
                     {"target": key, "decoys_per_active": args.decoys_per_active,
                      "max_decoy_similarity": MAX_DECOY_SIMILARITY})
    logger.info(f"Report written: {out_json.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
