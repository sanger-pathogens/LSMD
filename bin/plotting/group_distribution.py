#!/usr/bin/env python3
"""
Group/lineage composition chart -- sets up the "why specificity filtering matters"
story: here's the target group, and here's everything else a real marker has to
NOT cross-react with.

Design (dataviz method): horizontal bar, log-x (range here is target~5,000 down to
several groups at 3 -- linear would make everything but the target look like zero),
sorted by count descending. Only two colors, not one per group: blue for the target
group, neutral grey for every other group -- the point is "target vs. everything
else," not "22 individually distinct things," and a hue per group would both violate
the don't-cycle-past-~8-categorical-colors rule and dilute the actual story.

Usage: python3 group_distribution.py <metadata.csv> <target_group> [-o OUTDIR] [--group-col lineage]
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE = "#2a78d6"
GREY = "#9a9a96"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.size": 11,
        "font.family": "DejaVu Sans",
        "axes.edgecolor": INK2,
        "axes.linewidth": 0.8,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
    }
)


def main():
    p = argparse.ArgumentParser(description="Group/lineage composition bar chart")
    p.add_argument("metadata_csv")
    p.add_argument("target_group", help="the group to highlight, e.g. 7PET")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: current dir)")
    p.add_argument("--group-col", default="lineage", help="Metadata column holding the group label (default: lineage)")
    p.add_argument("--species", default="V. cholerae", help="Species name for the title (default: 'V. cholerae')")
    args = p.parse_args()

    with open(args.metadata_csv) as f:
        rows = list(csv.DictReader(f))
    counts = Counter(r[args.group_col].strip() for r in rows if r.get(args.group_col, "").strip())
    if args.target_group not in counts:
        sys.exit(f"target group '{args.target_group}' not found in {args.group_col} column")

    ordered = sorted(counts.items(), key=lambda kv: -kv[1])
    labels = [k for k, _ in ordered]
    values = [v for _, v in ordered]
    colors = [BLUE if lab == args.target_group else GREY for lab in labels]

    y = list(range(len(labels)))[::-1]
    fig, ax = plt.subplots(figsize=(10, 0.32 * len(labels) + 2))
    for yi, val, color in zip(y, values, colors):
        ax.barh(yi, val, color=color, alpha=0.9, edgecolor=INK2, height=0.65)
        ax.text(val * 1.15, yi, f"{val:,}", va="center", ha="left", fontsize=9, color=INK)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("Genomes (log scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Group", fontsize=11, fontweight="bold")
    # TITLE TBC -- suggestion, meant to be improved:
    ax.set_title(
        f"{args.species} Dataset Composition: {args.target_group} vs. {len(labels)-1} Other Groups",
        fontsize=13,
        fontweight="bold",
    )
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    total = sum(values)
    fig.text(
        0.5,
        -0.01,
        f"{total:,} genomes total across {len(labels)} groups -- blue = target group for marker discovery, "
        "grey = every group markers must not cross-react with",
        ha="center",
        va="top",
        fontsize=8,
        color=INK2,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{args.species.replace(' ', '_').replace('.', '')}_{args.target_group}_group_distribution.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
