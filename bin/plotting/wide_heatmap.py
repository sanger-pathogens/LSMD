#!/usr/bin/env python3
"""
Marker specificity heatmap from a wide TSV (one column per colour/lineage, one
row per marker) -- the shape lineage-level checks come in, as opposed to
marker_heatmap.py's long validation.tsv/group_detail.tsv pair from the ATB check.

Same visual convention as marker_heatmap.py: markers on x, colours on y (target
column pinned top, bold, solid), sequential white->blue fill for magnitude, red
border only on off-target cells at/above --border-min-frac.

Usage: python3 wide_heatmap.py <wide.tsv> <target_column> [-o OUTDIR] [--out-name NAME]
           [--border-min-frac 5.0] [--title TITLE] [--caption TEXT]
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"
RED = "#c0392b"
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", ["#fcfcfb", "#2a78d6"])

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.size": 11,
        "font.family": "DejaVu Sans",
        "text.color": INK,
    }
)


def main():
    p = argparse.ArgumentParser(description="Marker specificity heatmap from a wide (one-column-per-colour) TSV")
    p.add_argument("wide_tsv", help="header row = colour names, one data row per marker")
    p.add_argument("target_column", help="the on-target colour, e.g. 7PET")
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--out-name", default=None, help="output filename (default: <target_column>_wide_heatmap.png)")
    p.add_argument("--border-min-frac", type=float, default=5.0)
    p.add_argument("--title", default=None)
    p.add_argument("--caption", default=None, help="extra line drawn under the figure")
    args = p.parse_args()

    with open(args.wide_tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        sys.exit(f"No rows in {args.wide_tsv}")
    columns_in_file = list(rows[0].keys())
    if args.target_column not in columns_in_file:
        sys.exit(f"{args.target_column!r} not a column in {args.wide_tsv} (have {columns_in_file})")

    other_cols = [c for c in columns_in_file if c != args.target_column]
    columns = [args.target_column] + other_cols
    marker_order = [str(i + 1) for i in range(len(rows))]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / (args.out_name or f"{args.target_column}_wide_heatmap.png")

    n_rows, n_cols = len(columns), len(marker_order)
    fig, ax = plt.subplots(figsize=(max(6, 0.12 * n_cols + 2), max(3, 0.8 * n_rows + 1.5)))

    for yi, col in enumerate(columns):
        for xi, r in enumerate(rows):
            val = float(r[col])
            ax.add_patch(Rectangle((xi, yi), 1, 1, facecolor=SEQ_CMAP(val / 100.0), edgecolor=GRID, linewidth=0.4))
            if col != args.target_column and val >= args.border_min_frac:
                ax.add_patch(Rectangle((xi, yi), 1, 1, facecolor="none", edgecolor=RED, linewidth=1.3))

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_xticks([])
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels(columns, fontsize=10)
    for tick, col in zip(ax.get_yticklabels(), columns):
        if col == args.target_column:
            tick.set_fontweight("bold")
    ax.set_xlabel(f"{args.target_column} Markers (n={n_cols})", fontsize=11, fontweight="bold")
    ax.set_title(args.title or f"Candidate Markers vs Lineages/Sister Species: {args.target_column} Specificity",
                 fontsize=13, fontweight="bold")
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=SEQ_CMAP, norm=plt.Normalize(0, 100))
    cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("% of the marker's k-mers found in that colour", fontsize=9)

    caption = args.caption or f"Red outline = off-target colour hit at ≥{args.border_min_frac:g}% (not {args.target_column})"
    fig.text(0.5, -0.04, caption, ha="center", va="top", fontsize=8, color=INK2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
