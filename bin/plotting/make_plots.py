#!/usr/bin/env python3
"""Sprint-review plots for the marker-validation run.

Reads the per-tier outputs (vc_7PET_<mode>/tier_a|c/...) from a validation run
directory and writes slide-ready PNGs to <root>/plots/.

Usage: python3 make_plots.py [root]   -- root defaults to the current directory,
so run it from inside the validation run directory, or pass the path explicitly.

Tier C here is the STRICT blastn check (two knobs only: pident >= 98, qcovhsp >= 80).
A marker PASSes only if it is carried by >= 1 7PET genome and ZERO non-7PET genomes.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
OUT = ROOT / "plots"
OUT.mkdir(exist_ok=True)

ALL_MODES = ["core", "relaxed", "catchall"]
MODE_LABEL = {"core": "core", "relaxed": "relaxed", "catchall": "catchall"}
MODE_SUB = {"core": "within ≥ 0.95", "relaxed": "within ≥ 0.5", "catchall": "within > 0"}

# validated categorical palette (dataviz reference instance, light mode)
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 12, "font.family": "DejaVu Sans",
    "axes.edgecolor": INK2, "axes.linewidth": 0.8,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
})

TIER_C = {m: ROOT / f"vc_7PET_{m}" / "tier_c" / f"vc_7PET_{m}_validation.tsv" for m in ALL_MODES}
TALLY = {m: ROOT / f"vc_7PET_{m}" / "tier_c" / f"vc_7PET_{m}_lineage_tally.tsv" for m in ALL_MODES}
MODES = [m for m in ALL_MODES if TIER_C[m].exists()]


def tsv(path):
    return pd.read_csv(path, sep="\t")


def tier_a(mode):
    return tsv(ROOT / f"vc_7PET_{mode}" / "tier_a" / f"vc_7PET_{mode}_validation.tsv")


def header(fig, title, subtitle, y=0.965):
    fig.text(0.015, y, title, ha="left", va="top", fontsize=14.5, fontweight="bold")
    fig.text(0.015, y - 0.058, subtitle, ha="left", va="top", fontsize=9.5, color=INK2)


def style_ax(ax, axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)
    getattr(ax, f"{axis}axis").grid(True, color=GRID, linewidth=0.8)
    getattr(ax, f"{'x' if axis == 'y' else 'y'}axis").grid(False)


# ============================================================================
# Figure 1 - outcome of the strict blastn lineage check, by candidate mode
# ============================================================================
GREEN = "#008300"
GREY = "#a9a89f"


def fig_pass_rate():
    cats = []   # (flag%, rare%, usable%) per mode
    n = []
    for m in MODES:
        d = tsv(TIER_C[m])
        nm = len(d)
        n.append(nm)
        flag = (d["verdict"] == "FLAG").mean() * 100
        p = d["verdict"] == "PASS"
        usable = (p & (d["within_7PET"] >= 50)).mean() * 100
        rare = (p & (d["within_7PET"] < 50)).mean() * 100
        cats.append((flag, rare, usable))

    x = np.arange(len(MODES))
    w = 0.5
    fig = plt.figure(figsize=(9, 5.6))
    ax = fig.add_axes([0.10, 0.17, 0.87, 0.56])
    flags = [c[0] for c in cats]
    rares = [c[1] for c in cats]
    usables = [c[2] for c in cats]
    ax.bar(x, flags, w, color=ORANGE, edgecolor=SURFACE, linewidth=2,
           label="also carried by ≥ 1 non-7PET genome")
    ax.bar(x, rares, w, bottom=flags, color=GREY, edgecolor=SURFACE, linewidth=2,
           label="7PET-only, but in < 50% of 7PET  (too rare to be a marker)")
    ax.bar(x, usables, w, bottom=[f + r for f, r in zip(flags, rares)], color=GREEN,
           edgecolor=SURFACE, linewidth=2, label="7PET-only AND in ≥ 50% of 7PET  (candidate marker)")
    for xi, (f, r, u) in zip(x, cats):
        if f > 6:
            ax.annotate(f"{f:.0f}%", (xi, f / 2), ha="center", va="center", fontsize=10.5,
                        fontweight="bold", color="white")
        if r > 6:
            ax.annotate(f"{r:.0f}%", (xi, f + r / 2), ha="center", va="center", fontsize=10.5,
                        fontweight="bold", color="white")
        ax.annotate(f"{u:.1f}%", (xi, 103), ha="center", fontsize=10, fontweight="bold", color=GREEN)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{MODE_LABEL[m]}\n({MODE_SUB[m]})\nn = {n[i]:,}" for i, m in enumerate(MODES)])
    ax.set_ylabel("candidate markers  (%)")
    ax.set_ylim(0, 108)
    ax.set_yticks([0, 25, 50, 75, 100])
    style_ax(ax)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02), fontsize=9,
              handlelength=1.3, labelspacing=0.35)
    header(fig, "Not one candidate marker is both 7PET-private and 7PET-prevalent",
           "blastn vs the run's own 6,147-genome panel (pident ≥ 98, cov ≥ 80). The exact-31-mer "
           "specificity filter passed 100% of every set.")
    fig.savefig(OUT / "01_pass_rate_by_tier.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / "01_pass_rate_by_tier.png")


# ============================================================================
# Figure 2 - top non-7PET lineages the markers land in, sorted
# ============================================================================
def fig_lineages_matched():
    panes = [m for m in ALL_MODES if TALLY[m].exists()]
    if not panes:
        return
    fig, axes = plt.subplots(1, len(panes), figsize=(5.8 * len(panes), 6.5))
    if len(panes) == 1:
        axes = [axes]
    fig.subplots_adjust(left=0.09, right=0.98, top=0.70, bottom=0.13, wspace=0.62)

    TOP = 12
    for ax, m in zip(axes, panes):
        t = tsv(TALLY[m]).head(TOP).iloc[::-1]        # most-hit at the top
        y = np.arange(len(t))
        ax.barh(y, t["pct_markers_matched"], 0.66, color=ORANGE, edgecolor=SURFACE, linewidth=1.2,
                label="% of candidate markers hitting the lineage", zorder=2)
        ax.scatter(t["mean_frac_of_group"], y, s=34, color=BLUE, zorder=3,
                   label="mean % of that lineage's genomes hit")
        for yi, pct in enumerate(t["pct_markers_matched"]):
            ax.annotate(f"{pct:.0f}%", (pct, yi), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=9, fontweight="bold", color=INK)
        ax.set_yticks(y)
        ax.set_yticklabels(t["group"], fontsize=9.5)
        ax.set_xlim(0, 108)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("percent")
        n = len(tsv(TIER_C[m]))
        ax.set_title(f"{MODE_LABEL[m]}   (n = {n:,})", fontsize=12, fontweight="bold", pad=8)
        style_ax(ax, axis="x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    h, lab = axes[0].get_legend_handles_labels()
    order = [lab.index("% of candidate markers hitting the lineage"),
             lab.index("mean % of that lineage's genomes hit")]
    axes[0].legend([h[i] for i in order], [lab[i] for i in order], frameon=False,
                   loc="lower left", bbox_to_anchor=(0.0, 1.10), fontsize=9, handletextpad=0.5)
    header(fig, "Which other V. cholerae lineages carry the “7PET-specific” markers",
           "blastn pident ≥ 98 / coverage ≥ 80, sorted by prevalence. core/relaxed: every candidate also hits ≥ 1 "
           "sibling lineage, often across the whole lineage. catchall: ~60% do; the rest are just rare.", y=0.955)
    fig.savefig(OUT / "02_lineages_matched.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / "02_lineages_matched.png")


# ============================================================================
# Figure 3 - how much of the species carries each marker (strict blastn)
# ============================================================================
def fig_genome_spread():
    means = [tsv(TIER_C[m])["n_genomes_hit"].mean() for m in MODES]
    frac = [100 * v / 6147 for v in means]

    fig = plt.figure(figsize=(8.6, 5.0))
    ax = fig.add_axes([0.16, 0.16, 0.80, 0.58])
    xs = np.arange(len(MODES))
    bars = ax.bar(xs, frac, 0.5, color=BLUE, edgecolor=SURFACE, linewidth=1.5)
    for r, mv in zip(bars, means):
        ax.annotate(f"{r.get_height():.0f}%   ({mv:,.0f} / 6,147)",
                    (r.get_x() + r.get_width() / 2, r.get_height()),
                    xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                    fontsize=10.5, fontweight="bold")
    ax.axhline(100 * 5001 / 6147, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)))
    ax.annotate("7PET is ≈ 81% of the panel", (xs[-1] + 0.35, 100 * 5001 / 6147),
                xytext=(0, -13), textcoords="offset points", ha="right", fontsize=8.5, color=INK2)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{MODE_LABEL[m]}\n({MODE_SUB[m]})" for m in MODES])
    ax.set_xlim(-0.6, len(MODES) - 0.4)
    ax.set_ylabel("mean V. cholerae assemblies\ncarrying each marker  (%)")
    ax.set_ylim(0, 118)
    ax.set_yticks([0, 25, 50, 75, 100])
    style_ax(ax)
    header(fig, "How much of the species each candidate marker covers",
           "Mean genomes (of 6,147) with a near-identical, near-full-length blastn hit "
           "(pident ≥ 98, cov ≥ 80). core sits above 7PET's own\nsize — it spans sibling lineages too; "
           "catchall candidates are mostly too rare to matter.")
    fig.savefig(OUT / "03_genomes_carrying_marker.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / "03_genomes_carrying_marker.png")


if __name__ == "__main__":
    print("modes with Tier C output:", ", ".join(MODES) or "(none)")
    fig_pass_rate()
    fig_lineages_matched()
    fig_genome_spread()
