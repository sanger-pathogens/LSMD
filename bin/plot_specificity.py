#!/usr/bin/env python3
"""
Specificity score figures for a core_catchall_filter.py E output.

score = core_pct - outside_pct (Youden's J): +1 = perfect marker (100% of the
target lineage, 0% of everything else); 0 = no discriminating power; negative
= the unitig is actually more common outside the lineage than inside it.

Usage:
    python3 plot_specificity.py --specificity <lineage>_<mode>_specificity.tsv \
        --label <e.g. v_cholerae_sublineage_1.0_E> --display-name <e.g. 'V. cholerae -- lineage 1.0'> \
        --mode core --out-dir <dir>
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# Palette (dataviz skill reference palette, series slot 1 -- light mode)
SERIES_BLUE = "#2a78d6"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e3e2dd"
SURFACE = "#fcfcfb"
THRESHOLD_COLOR = "#e34948"  # palette slot 8, red -- reserved/status-style reference line

GOOD_SPECIFICITY_THRESHOLD = 0.9

PREVIEW_NOTE = (
    "PREVIEW on E (lineage-core candidates) — cross-lineage exclusion "
    "(F = E − D) not yet applied"
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--specificity", required=True, type=Path,
                   help="'{lineage_id}_{mode}_specificity.tsv' from core_catchall_filter.py "
                        "(columns: unitig_id, core_pct, outside_pct, specificity_score)")
    p.add_argument("--label", required=True,
                   help="Short identifier for filenames, e.g. v_cholerae_sublineage_1.0_E")
    p.add_argument("--display-name", required=True,
                   help="Human-readable species/lineage name for chart titles, "
                        "e.g. 'V. cholerae — lineage 1.0'")
    p.add_argument("--mode", required=True, choices=["core", "relaxed", "catchall"],
                   help="Threshold mode used to build E, shown on every figure so "
                        "it's unambiguous if the figure is shared without context")
    p.add_argument("--good-threshold", type=float, default=GOOD_SPECIFICITY_THRESHOLD,
                   help=f"Specificity score marking a 'high-specificity' candidate marker "
                        f"(default {GOOD_SPECIFICITY_THRESHOLD})")
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=TEXT_SECONDARY, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(TEXT_PRIMARY)
    ax.yaxis.label.set_color(TEXT_PRIMARY)


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    threshold = args.good_threshold

    data = np.loadtxt(args.specificity, skiprows=1, usecols=(1, 2, 3))
    core_pct, outside_pct, score = data[:, 0], data[:, 1], data[:, 2]
    n = len(score)
    n_good = int((score >= threshold).sum())
    n_neg = int((score < 0).sum())

    stats_lines = [
        f"n unitigs (E, {args.mode} mode) : {n:,}",
        f"core_pct range          : {core_pct.min():.3f} - {core_pct.max():.3f}",
        f"outside_pct range       : {outside_pct.min():.3f} - {outside_pct.max():.3f}",
        f"specificity_score min   : {score.min():.4f}",
        f"specificity_score median: {np.median(score):.4f}",
        f"specificity_score mean  : {score.mean():.4f}",
        f"specificity_score max   : {score.max():.4f}",
        f"score >= {threshold:g} (high-specificity)      : {n_good:,} / {n:,} ({100 * n_good / n:.2f}%)",
        f"score < 0 (more common outside than in)  : {n_neg:,} / {n:,} ({100 * n_neg / n:.2f}%)",
    ]
    print("\n".join(stats_lines))
    (args.out_dir / f"{args.label}_specificity_stats.txt").write_text(
        "\n".join(stats_lines) + "\n"
    )

    fig_fmt = dict(figsize=(9, 4.6), dpi=200)

    # --- Figure 1: log-count histogram, full range -----------------------------
    # Almost all candidates sit in a narrow spike around score=0 -- they pass
    # the "core" (within-lineage) threshold but are just as common outside the
    # lineage, so they carry no discriminating power. The unitigs that actually
    # make good diagnostic markers are a thin tail out near score=1. A
    # linear-count histogram renders that tail as invisible; log-count keeps
    # the near-zero spike AND resolves the tail in the same panel, without
    # touching the x-axis (score is bounded, not heavy-tailed itself -- only
    # the counts are).
    fig, ax = plt.subplots(**fig_fmt)
    bins = np.linspace(score.min(), score.max(), 80)
    ax.hist(score, bins=bins, color=SERIES_BLUE, edgecolor=SURFACE, linewidth=0.4, zorder=2)
    ax.set_yscale("log")
    ax.axvline(threshold, color=THRESHOLD_COLOR, linewidth=1.5,
               linestyle=(0, (4, 2)), zorder=3)
    ax.text(threshold, ax.get_ylim()[1], f" score ≥ {threshold:g} = high-specificity candidate",
            color=THRESHOLD_COLOR, fontsize=9, va="top", ha="left")
    ax.set_xlabel("Specificity score (core_pct − outside_pct)")
    ax.set_ylabel("Count (log scale)")
    ax.set_title(
        f"{args.display_name} — candidate marker specificity, {args.mode} mode (n={n:,})",
        color=TEXT_PRIMARY, fontsize=11, loc="left")
    ax.text(0.0, 1.10, PREVIEW_NOTE, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, color=TEXT_SECONDARY, style="italic")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(args.out_dir / f"{args.label}_specificity_score_hist.png", facecolor=SURFACE)
    plt.close(fig)

    # --- Figure 2: survival curve (fraction clearing a given score) ------------
    # Answers the practical question directly: what fraction of E's candidates
    # are actually usable as lineage-specific markers, at any given cutoff.
    sorted_scores = np.sort(score)
    cdf = np.arange(1, n + 1) / n
    survival = 1.0 - cdf
    frac_good = n_good / n

    fig, ax = plt.subplots(**fig_fmt)
    ax.plot(sorted_scores, survival, color=SERIES_BLUE, linewidth=2, zorder=2)
    ax.axvline(threshold, color=THRESHOLD_COLOR, linewidth=1.5,
               linestyle=(0, (4, 2)), zorder=3)
    ax.axhline(frac_good, color=THRESHOLD_COLOR, linewidth=1, alpha=0.5,
               linestyle=(0, (1, 2)), zorder=1)
    ax.annotate(f"{frac_good * 100:.2f}% of candidates are high-specificity (score ≥ {threshold:g})",
                xy=(threshold, frac_good), xytext=(-10, 25), ha="right",
                textcoords="offset points", fontsize=9, color=THRESHOLD_COLOR,
                arrowprops=dict(arrowstyle="-", color=THRESHOLD_COLOR, lw=0.8))
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))
    ax.set_xlabel("Specificity score cutoff (core_pct − outside_pct)")
    ax.set_ylabel("Fraction of candidates ≥ cutoff")
    ax.set_title(
        f"{args.display_name} — candidate markers clearing a specificity cutoff, {args.mode} mode (n={n:,})",
        color=TEXT_PRIMARY, fontsize=11, loc="left")
    ax.text(0.0, 1.10, PREVIEW_NOTE, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, color=TEXT_SECONDARY, style="italic")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(args.out_dir / f"{args.label}_specificity_survival.png", facecolor=SURFACE)
    plt.close(fig)

    print(f"\nWrote figures + stats to {args.out_dir}")


if __name__ == "__main__":
    main()
