#!/usr/bin/env python3
"""
Specificity score figures for a lineage_specificity_score.py output.

Input is '{lineage_id}_specificity.tsv' -- one row per unitig in the
SPECIES-WIDE export graph (NOT the lineage-core candidate set E from
core_catchall_filter.py, which is an independently-built graph with different
unitig IDs). It's a self-consistent species-wide diagnostic view of how
discriminating each unitig would be for this lineage.

score = within_pct - outside_pct (Youden's J): +1 = perfect marker (100% of the
target lineage, 0% of everything else); 0 = no discriminating power; negative
= the unitig is actually more common outside the lineage than inside it.

Usage:
    python3 plot_specificity.py --specificity <lineage_id>_specificity.tsv \
        --label <e.g. v_cholerae_sublineage_1.0> --display-name <e.g. 'V. cholerae -- lineage 1.0'> \
        --mode core --out-dir <dir>
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

# Palette (dataviz skill reference palette, series slot 1 -- light mode)
SERIES_BLUE = "#2a78d6"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e3e2dd"
SURFACE = "#fcfcfb"
THRESHOLD_COLOR = "#e34948"  # palette slot 8, red -- reserved/status-style reference line

GOOD_SPECIFICITY_THRESHOLD = 0.9

PREVIEW_NOTE = (
    "Species-wide diagnostic — every unitig in the species graph scored " "for this lineage; not the candidate set E"
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--specificity",
        required=True,
        type=Path,
        help="'{lineage_id}_specificity.tsv' from lineage_specificity_score.py "
        "(columns: unitig_id, within_pct, outside_pct, specificity_score); "
        "one row per species-wide export unitig",
    )
    p.add_argument("--label", required=True, help="Short identifier for filenames, e.g. v_cholerae_sublineage_1.0")
    p.add_argument(
        "--display-name",
        required=True,
        help="Human-readable species/lineage name for chart titles, " "e.g. 'V. cholerae — lineage 1.0'",
    )
    p.add_argument(
        "--mode",
        required=True,
        choices=["core", "relaxed", "catchall"],
        help="within-lineage threshold mode of the companion core_catchall_filter "
        "run, shown on every figure for context if it's shared without it",
    )
    p.add_argument(
        "--good-threshold",
        type=float,
        default=GOOD_SPECIFICITY_THRESHOLD,
        help=f"Specificity score marking a 'high-specificity' candidate marker "
        f"(default {GOOD_SPECIFICITY_THRESHOLD})",
    )
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
    within_pct, outside_pct, score = data[:, 0], data[:, 1], data[:, 2]
    n = len(score)
    n_good = int((score >= threshold).sum())
    n_neg = int((score < 0).sum())

    stats_lines = [
        f"n unitigs scored (species-wide graph) : {n:,}",
        f"within_pct range        : {within_pct.min():.3f} - {within_pct.max():.3f}",
        f"outside_pct range       : {outside_pct.min():.3f} - {outside_pct.max():.3f}",
        f"specificity_score min   : {score.min():.4f}",
        f"specificity_score median: {np.median(score):.4f}",
        f"specificity_score mean  : {score.mean():.4f}",
        f"specificity_score max   : {score.max():.4f}",
        f"score >= {threshold:g} (high-specificity)      : {n_good:,} / {n:,} ({100 * n_good / n:.2f}%)",
        f"score < 0 (more common outside than in)  : {n_neg:,} / {n:,} ({100 * n_neg / n:.2f}%)",
    ]
    print("\n".join(stats_lines))
    (args.out_dir / f"{args.label}_specificity_stats.txt").write_text("\n".join(stats_lines) + "\n")

    fig_fmt = dict(figsize=(9, 4.6), dpi=200)

    # --- Figure 1: linear-count histogram + zoomed inset on the marker tail ----
    # Almost all species-wide unitigs sit in a narrow spike around score=0 --
    # present within the lineage but just as common outside it, so they carry
    # no discriminating power. The unitigs that make good diagnostic markers are
    # a thin tail out near score=1. A log y-axis resolves that tail but badly
    # distorts how overwhelming the near-zero spike is (a reader eyeballs bar
    # heights linearly). Instead: honest linear-count main panel -- the spike
    # reads at true scale -- plus an inset zoomed to score >= good-threshold with
    # its own y-scale, so the handful of real candidates are still visible and
    # countable.
    fig, ax = plt.subplots(**fig_fmt)
    bins = np.linspace(score.min(), score.max(), 80)
    ax.hist(score, bins=bins, color=SERIES_BLUE, edgecolor=SURFACE, linewidth=0.4, zorder=2)
    ax.axvline(threshold, color=THRESHOLD_COLOR, linewidth=1.5, linestyle=(0, (4, 2)), zorder=3)
    ax.text(
        threshold - 0.015,
        ax.get_ylim()[1],
        f"score ≥ {threshold:g} = high-specificity ",
        color=THRESHOLD_COLOR,
        fontsize=8.5,
        va="top",
        ha="right",
        rotation=90,
    )
    ax.set_xlabel("Specificity score (within_pct − outside_pct)")
    ax.set_ylabel("Count")
    ax.set_title(
        f"{args.display_name} — species-wide unitig specificity, {args.mode} mode (n={n:,})",
        color=TEXT_PRIMARY,
        fontsize=11,
        loc="left",
    )
    ax.text(
        0.0,
        1.10,
        PREVIEW_NOTE,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=TEXT_SECONDARY,
        style="italic",
    )
    style_axes(ax)

    # Inset: the high-specificity tail only. Skipped if nothing clears the
    # threshold (an empty Axes reads as a bug, not as "zero markers").
    if n_good > 0:
        axin = ax.inset_axes([0.10, 0.42, 0.40, 0.46])
        tail = score[score >= threshold]
        axin.hist(
            tail,
            bins=np.linspace(threshold, score.max(), 20),
            color=THRESHOLD_COLOR,
            edgecolor=SURFACE,
            linewidth=0.4,
        )
        axin.set_title(f"score ≥ {threshold:g}  (n={n_good:,})", fontsize=8, color=TEXT_SECONDARY)
        axin.tick_params(labelsize=7, colors=TEXT_SECONDARY)
        for spine in ("top", "right"):
            axin.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.out_dir / f"{args.label}_specificity_score_hist.png", facecolor=SURFACE)
    plt.close(fig)

    # --- Figure 2: survival curve (fraction clearing a given score) ------------
    # Answers the practical question directly: what fraction of the species-wide
    # unitigs are usable as lineage-specific markers at any given cutoff --
    # i.e. marker yield vs specificity stringency.
    sorted_scores = np.sort(score)
    cdf = np.arange(1, n + 1) / n
    survival = 1.0 - cdf
    frac_good = n_good / n

    fig, ax = plt.subplots(**fig_fmt)
    ax.plot(sorted_scores, survival, color=SERIES_BLUE, linewidth=2, zorder=2)
    ax.axvline(threshold, color=THRESHOLD_COLOR, linewidth=1.5, linestyle=(0, (4, 2)), zorder=3)
    ax.axhline(frac_good, color=THRESHOLD_COLOR, linewidth=1, alpha=0.5, linestyle=(0, (1, 2)), zorder=1)
    ax.annotate(
        f"{frac_good * 100:.2f}% of unitigs are high-specificity (score ≥ {threshold:g})",
        xy=(threshold, frac_good),
        xytext=(-10, 25),
        ha="right",
        textcoords="offset points",
        fontsize=9,
        color=THRESHOLD_COLOR,
        arrowprops=dict(arrowstyle="-", color=THRESHOLD_COLOR, lw=0.8),
    )
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))
    ax.set_xlabel("Specificity score cutoff (within_pct − outside_pct)")
    ax.set_ylabel("Fraction of unitigs ≥ cutoff")
    ax.set_title(
        f"{args.display_name} — marker yield vs specificity stringency, {args.mode} mode (n={n:,})",
        color=TEXT_PRIMARY,
        fontsize=11,
        loc="left",
    )
    ax.text(
        0.0,
        1.10,
        PREVIEW_NOTE,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color=TEXT_SECONDARY,
        style="italic",
    )
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(args.out_dir / f"{args.label}_specificity_survival.png", facecolor=SURFACE)
    plt.close(fig)

    print(f"\nWrote figures + stats to {args.out_dir}")


if __name__ == "__main__":
    main()
