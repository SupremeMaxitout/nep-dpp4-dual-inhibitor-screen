#!/usr/bin/env python3
"""Generate every figure used in the project report."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

OUT = Path("figs"); OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

INK = "#1a1a1a"
NEP = "#2C5F7C"      # blue
DPP = "#A85751"      # red
NEUTRAL = "#6B6B6B"
ACCENT = "#3E7A5E"   # green
WARN = "#B07A2E"


# ===========================================================================
# Figure 1: pipeline overview flowchart
# ===========================================================================
def fig_pipeline():
    fig, ax = plt.subplots(figsize=(7.2, 8.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 24); ax.axis("off")

    stages = [
        (21.6, "Stage 1  Target preparation and validation",
         "Clean structures, retain catalytic zinc, derive grid box.\n"
         "Redock reference ligand. Benchmark against matched decoys.", ACCENT),
        (18.2, "Stage 2  Library construction",
         "Mine ChEMBL for measured binders. Filter on potency,\n"
         "assay confidence and physicochemical properties.", ACCENT),
        (14.8, "Stage 3  Primary docking and rescoring",
         "Dock each library against its own target. Rank by affinity\n"
         "and by ligand efficiency. Apply calibrated threshold.", ACCENT),
        (11.4, "Stage 4  Cross-docking and dual selection",
         "Dock each surviving set against the opposite target.\n"
         "Dual actives are the intersection.", ACCENT),
        (8.0, "Stage 5  Scaffold analysis and enumeration",
         "Extract shared scaffolds. Enforce dual pharmacophore.\n"
         "Enumerate analogues. Filter on synthetic accessibility.", ACCENT),
        (4.6, "Stage 6  Antitarget counter-screen",
         "Screen against ACE and DPP-8/9. Penalise cross-reactivity.", WARN),
        (1.9, "Prioritised set for in vitro testing",
         "Selectivity-filtered, with full provenance per compound.", NEP),
    ]

    for y, title, body, colour in stages:
        h = 2.5 if body else 1.6
        box = FancyBboxPatch((0.4, y - 0.1), 9.2, h,
                             boxstyle="round,pad=0.12,rounding_size=0.12",
                             linewidth=1.1, edgecolor=colour,
                             facecolor=colour, alpha=0.08)
        ax.add_patch(box)
        ax.text(0.75, y + h - 0.55, title, fontsize=9.5, fontweight="bold",
                color=colour, va="center")
        if body:
            ax.text(0.75, y + 0.62, body, fontsize=7.8, color=INK,
                    va="center", linespacing=1.45)

    arrows = [(21.5, 20.85), (18.1, 17.45), (14.7, 14.05), (11.3, 10.65),
              (7.9, 7.25), (4.45, 4.15)]
    for top, bottom in arrows:
        ax.annotate("", xy=(5, bottom), xytext=(5, top),
                    arrowprops=dict(arrowstyle="-|>", linewidth=1.3,
                                    color=NEUTRAL, shrinkA=2, shrinkB=2))

    plt.savefig(OUT / "fig1_pipeline.pdf")
    plt.close()


# ===========================================================================
# Figure 2: the pharmacophore conflict
# ===========================================================================
def fig_pharmacophore():
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.set_xlim(0, 12); ax.set_ylim(0, 6); ax.axis("off")

    # NEP side
    ax.add_patch(FancyBboxPatch((0.3, 2.6), 4.4, 3.0,
                 boxstyle="round,pad=0.15,rounding_size=0.15",
                 linewidth=1.3, edgecolor=NEP, facecolor=NEP, alpha=0.08))
    ax.text(2.5, 5.15, "Neprilysin", fontsize=11, fontweight="bold",
            color=NEP, ha="center")
    ax.text(2.5, 4.62, "zinc metallopeptidase", fontsize=7.5,
            color=NEUTRAL, ha="center", style="italic")
    ax.text(2.5, 3.95, "Catalytic Zn(II) held by\nHis583, His587, Glu646",
            fontsize=8, color=INK, ha="center", linespacing=1.4)
    ax.text(2.5, 3.05, "Requires an ANIONIC group\n(carboxylate, thiol, hydroxamate)",
            fontsize=8.2, color=NEP, ha="center", fontweight="bold",
            linespacing=1.4)

    # DPP side
    ax.add_patch(FancyBboxPatch((7.3, 2.6), 4.4, 3.0,
                 boxstyle="round,pad=0.15,rounding_size=0.15",
                 linewidth=1.3, edgecolor=DPP, facecolor=DPP, alpha=0.08))
    ax.text(9.5, 5.15, "DPP-IV", fontsize=11, fontweight="bold",
            color=DPP, ha="center")
    ax.text(9.5, 4.62, "serine protease", fontsize=7.5,
            color=NEUTRAL, ha="center", style="italic")
    ax.text(9.5, 3.95, "Glu205 / Glu206 dyad\nin the S2 subsite",
            fontsize=8, color=INK, ha="center", linespacing=1.4)
    ax.text(9.5, 3.05, "Requires a BASIC nitrogen\n(primary or secondary amine)",
            fontsize=8.2, color=DPP, ha="center", fontweight="bold",
            linespacing=1.4)

    # centre
    ax.add_patch(FancyBboxPatch((4.95, 3.35), 2.1, 1.6,
                 boxstyle="round,pad=0.12,rounding_size=0.12",
                 linewidth=1.4, edgecolor=WARN, facecolor="white"))
    ax.text(6.0, 4.42, "One molecule", fontsize=8.5, ha="center",
            fontweight="bold", color=INK)
    ax.text(6.0, 3.85, "must be acidic\nAND basic", fontsize=8,
            ha="center", color=WARN, linespacing=1.4)

    for x0, x1 in [(4.75, 4.95), (7.05, 7.25)]:
        ax.annotate("", xy=(x1, 4.15), xytext=(x0, 4.15),
                    arrowprops=dict(arrowstyle="-|>", linewidth=1.4,
                                    color=NEUTRAL))

    ax.text(6.0, 1.75, "Consequence: the molecule is a zwitterion at physiological pH.",
            fontsize=8.6, ha="center", color=INK, fontweight="bold")
    ax.text(6.0, 1.15,
            "High polarity and low membrane permeability follow directly.\n"
            "Sacubitril is administered as an ester prodrug for exactly this reason.",
            fontsize=8, ha="center", color=NEUTRAL, linespacing=1.5)

    plt.savefig(OUT / "fig2_pharmacophore.pdf")
    plt.close()


# ===========================================================================
# Figure 3: ROC curves, illustrative
# ===========================================================================
def fig_roc():
    from sklearn.metrics import roc_curve, roc_auc_score
    rng = np.random.default_rng(42)
    n, n_act = 1000, 50
    labels = np.zeros(n, int); labels[:n_act] = 1; rng.shuffle(labels)

    scenarios = [
        ("Strong protocol", 1.8, ACCENT, "-"),
        ("Acceptable protocol", 1.2, NEP, "-"),
        ("Weak protocol", 0.5, WARN, "--"),
        ("Random", 0.0, NEUTRAL, ":"),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2))

    for label, sep, colour, ls in scenarios:
        scores = np.where(labels == 1, -8.0 - sep, -8.0) + rng.normal(0, 1.2, n)
        fpr, tpr, _ = roc_curve(labels, -scores)
        auc = roc_auc_score(labels, -scores)
        ax1.plot(fpr, tpr, color=colour, linestyle=ls, linewidth=1.6,
                 label=f"{label}  (AUC {auc:.2f})")

    ax1.plot([0, 1], [0, 1], color="#cccccc", linewidth=0.8, zorder=0)
    ax1.set_xlabel("False positive rate", fontsize=8.5)
    ax1.set_ylabel("True positive rate", fontsize=8.5)
    ax1.set_title("ROC behaviour across protocol quality", fontsize=9,
                  fontweight="bold", pad=8)
    ax1.legend(fontsize=6.8, loc="lower right", frameon=False)
    ax1.tick_params(labelsize=7.5)
    ax1.spines[["top", "right"]].set_visible(False)

    # EF bar chart
    def ef(labels, scores, frac):
        n = len(labels); ntop = max(1, int(round(n * frac)))
        order = np.argsort(scores)
        return (labels[order][:ntop].sum() / ntop) / (labels.sum() / n)

    fracs = [0.005, 0.01, 0.05]
    width = 0.2
    xs = np.arange(len(fracs))
    for i, (label, sep, colour, _) in enumerate(scenarios):
        scores = np.where(labels == 1, -8.0 - sep, -8.0) + rng.normal(0, 1.2, n)
        vals = [ef(labels, scores, f) for f in fracs]
        ax2.bar(xs + i * width - 1.5 * width, vals, width,
                color=colour, alpha=0.85, label=label)

    ax2.axhline(1.0, color="#999999", linewidth=0.8, linestyle="--")
    ax2.text(2.35, 1.35, "chance", fontsize=6.5, color="#999999", ha="right")
    ax2.set_xticks(xs)
    ax2.set_xticklabels(["EF 0.5%", "EF 1%", "EF 5%"], fontsize=8)
    ax2.set_ylabel("Enrichment factor", fontsize=8.5)
    ax2.set_title("Early enrichment", fontsize=9, fontweight="bold", pad=8)
    ax2.tick_params(labelsize=7.5)
    ax2.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(OUT / "fig3_roc.pdf")
    plt.close()


# ===========================================================================
# Figure 4: attrition funnel
# ===========================================================================
def fig_funnel():
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.set_xlim(0, 10); ax.set_ylim(0, 9); ax.axis("off")

    stages = [
        ("ChEMBL bioactivity records retrieved", 100.0, "~250,000", NEUTRAL),
        ("Unique compounds after aggregation", 62.0, "~40,000", NEUTRAL),
        ("Survive cleaning and property filters", 44.0, "~22,000", NEP),
        ("Pass primary docking threshold", 21.0, "~1,800", NEP),
        ("Active against BOTH targets (cross-dock)", 9.5, "~120", ACCENT),
        ("Carry a valid dual pharmacophore", 5.5, "~35", ACCENT),
        ("Survive ACE / DPP-8/9 counter-screen", 3.0, "~12", WARN),
    ]

    y = 8.3
    for name, width_pct, count, colour in stages:
        w = width_pct / 100 * 7.2
        x0 = 1.1 + (7.2 - w) / 2
        ax.add_patch(FancyBboxPatch((x0, y - 0.42), w, 0.68,
                     boxstyle="round,pad=0.03,rounding_size=0.06",
                     linewidth=0.9, edgecolor=colour, facecolor=colour,
                     alpha=0.22))
        ax.text(x0 + w / 2, y - 0.08, count, fontsize=8, ha="center",
                va="center", fontweight="bold", color=INK)
        ax.text(1.0, y - 0.08, name, fontsize=7.6, ha="right", va="center",
                color=INK)
        y -= 1.15

    ax.text(5.0, 0.35,
            "Order-of-magnitude attrition for a screen of this design. Absolute numbers depend on\n"
            "target, thresholds and library composition, and are illustrative rather than measured.",
            fontsize=7, ha="center", color=NEUTRAL, style="italic", linespacing=1.5)

    plt.savefig(OUT / "fig4_funnel.pdf")
    plt.close()


# ===========================================================================
# Figure 5: ML concept, docking-score label vs experimental label
# ===========================================================================
def fig_ml():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    rng = np.random.default_rng(7)

    # Panel A: docking vs experiment, weak correlation
    n = 180
    true_p = rng.normal(7.0, 1.2, n)
    dock = -(4.4 + 0.42 * true_p) + rng.normal(0, 0.95, n)
    ax = axes[0]
    ax.scatter(true_p, dock, s=9, color=NEP, alpha=0.55, edgecolors="none")
    r = np.corrcoef(true_p, dock)[0, 1]
    ax.set_xlabel("Experimental pChEMBL", fontsize=8)
    ax.set_ylabel("Docking score (kcal/mol)", fontsize=8)
    ax.set_title(f"Docking vs experiment\nr = {r:.2f}", fontsize=8.5,
                 fontweight="bold", pad=6)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    # Panel B: what a model trained on docking scores learns
    ax = axes[1]
    pred_dock = dock + rng.normal(0, 0.28, n)
    ax.scatter(dock, pred_dock, s=9, color=WARN, alpha=0.55, edgecolors="none")
    lim = [dock.min() - 0.4, dock.max() + 0.4]
    ax.plot(lim, lim, color="#999999", linewidth=0.8, linestyle="--")
    r2 = np.corrcoef(dock, pred_dock)[0, 1] ** 2
    ax.set_xlabel("Docking score (label)", fontsize=8)
    ax.set_ylabel("Model prediction", fontsize=8)
    ax.set_title(f"Trained on docking scores\nR$^2$ = {r2:.2f} vs the WRONG target",
                 fontsize=8.5, fontweight="bold", pad=6)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    # Panel C: same model against experiment
    ax = axes[2]
    implied = (-(pred_dock) - 4.4) / 0.42
    ax.scatter(true_p, implied, s=9, color=DPP, alpha=0.55, edgecolors="none")
    lim = [true_p.min() - 0.3, true_p.max() + 0.3]
    ax.plot(lim, lim, color="#999999", linewidth=0.8, linestyle="--")
    r2b = np.corrcoef(true_p, implied)[0, 1] ** 2
    ax.set_xlabel("Experimental pChEMBL", fontsize=8)
    ax.set_ylabel("Implied prediction", fontsize=8)
    ax.set_title(f"Same model vs reality\nR$^2$ = {r2b:.2f}", fontsize=8.5,
                 fontweight="bold", pad=6)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(OUT / "fig5_ml.pdf")
    plt.close()


# ===========================================================================
# Figure 6: multi-task architecture
# ===========================================================================
def fig_multitask():
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    ax.set_xlim(0, 12); ax.set_ylim(0, 5.5); ax.axis("off")

    def box(x, y, w, h, text, colour, fs=7.8, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                     boxstyle="round,pad=0.08,rounding_size=0.1",
                     linewidth=1.1, edgecolor=colour, facecolor=colour,
                     alpha=0.12))
        ax.text(x + w / 2, y + h / 2, text, fontsize=fs, ha="center",
                va="center", color=INK, linespacing=1.4,
                fontweight="bold" if bold else "normal")

    box(0.2, 2.1, 1.9, 1.3, "Molecule\n(SMILES)", NEUTRAL)
    box(2.5, 2.1, 2.1, 1.3, "Featurisation\nECFP4 +\ndescriptors", NEUTRAL)
    box(5.0, 2.1, 2.2, 1.3, "Shared\nrepresentation", NEP, bold=True)
    box(7.9, 3.3, 2.0, 1.0, "Head 1\nNEP pKi", NEP)
    box(7.9, 1.2, 2.0, 1.0, "Head 2\nDPP-IV pKi", DPP)
    box(10.3, 2.2, 1.5, 1.1, "Dual\nscore", ACCENT, bold=True)

    conns = [(2.1, 2.75, 2.5, 2.75), (4.6, 2.75, 5.0, 2.75),
             (7.2, 2.9, 7.9, 3.8), (7.2, 2.6, 7.9, 1.7),
             (9.9, 3.8, 10.3, 3.0), (9.9, 1.7, 10.3, 2.5)]
    for x0, y0, x1, y1 in conns:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                     arrowstyle="-|>", mutation_scale=9,
                     linewidth=1.1, color=NEUTRAL))

    ax.text(6.0, 0.55,
            "The shared layer learns chemistry relevant to both targets. Each head specialises.\n"
            "Training on experimental pKi rather than docking scores keeps the labels tied to biology.",
            fontsize=7.6, ha="center", color=NEUTRAL, linespacing=1.6)

    plt.savefig(OUT / "fig6_multitask.pdf")
    plt.close()


# ===========================================================================
# Figure 7: RMSD symmetry correction
# ===========================================================================
def fig_rmsd():
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    methods = ["Naive\natom-index RMSD", "Symmetry-corrected\nGetBestRMS"]
    values = [0.995, 0.000]
    colours = [DPP, ACCENT]
    bars = ax.bar(methods, values, color=colours, alpha=0.85, width=0.55)
    ax.axhline(2.0, color=NEUTRAL, linestyle="--", linewidth=1.0)
    ax.text(1.45, 2.08, "2.0 A pass threshold", fontsize=6.8,
            color=NEUTRAL, ha="right")
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.08,
                f"{v:.3f} A", ha="center", fontsize=8, fontweight="bold")
    ax.set_ylabel("RMSD (Angstrom)", fontsize=8.5)
    ax.set_ylim(0, 2.5)
    ax.set_title("Identical pose, carboxylate oxygens exchanged",
                 fontsize=8.5, fontweight="bold", pad=8)
    ax.tick_params(labelsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(OUT / "fig7_rmsd.pdf")
    plt.close()


if __name__ == "__main__":
    fig_pipeline(); print("fig1 pipeline")
    fig_pharmacophore(); print("fig2 pharmacophore")
    fig_roc(); print("fig3 roc")
    fig_funnel(); print("fig4 funnel")
    fig_ml(); print("fig5 ml")
    fig_multitask(); print("fig6 multitask")
    fig_rmsd(); print("fig7 rmsd")
    print("\nAll figures written to figs/")
