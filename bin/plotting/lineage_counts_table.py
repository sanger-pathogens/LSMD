#!/usr/bin/env python3
"""
Lineage/group composition table -- the input-data-table companion to
group_distribution.py's bar chart, for slides where a table reads better than
a chart (e.g. an appendix/methods slide listing exact counts).

Usage: python3 lineage_counts_table.py <color_mapping_stats.json> [-o OUTDIR] [--species "V. cholerae"] [--highlight 7PET]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"
BLUE = "#2a78d6"
HEADER_BG = "#eef2f8"
HIGHLIGHT_BG = "#dce8fa"

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
    p = argparse.ArgumentParser(description="Lineage/group counts table from color_mapping's *_stats.json")
    p.add_argument("stats_json", help="<species>_stats.json from BUILD_COLOR_INDEX:COLOR_MAPPING")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: current dir)")
    p.add_argument("--species", default="V. cholerae", help="Species name for the title (default: 'V. cholerae')")
    p.add_argument("--highlight", default="7PET", help="Group to highlight as the marker-discovery target (default: 7PET)")
    p.add_argument("--exclude", default=None, help="Comma-separated group names to drop from the table entirely (e.g. 'unclassified' -- the no-metadata-row bucket, not a real lineage)")
    args = p.parse_args()

    with open(args.stats_json) as f:
        stats = json.load(f)
    groups = dict(stats["assemblies_per_group"])
    excluded = [g.strip() for g in args.exclude.split(",")] if args.exclude else []
    n_excluded = sum(groups.pop(g, 0) for g in excluded)
    total = sum(groups.values())
    ordered = sorted(groups.items(), key=lambda kv: -kv[1])

    n_rows = len(ordered) + 2  # header + rows + total
    fig, ax = plt.subplots(figsize=(6.5, 0.32 * n_rows + 0.3))
    ax.axis("off")

    col_labels = ["Group", "Genomes", "% of total"]
    cell_text = [[g, f"{n:,}", f"{100 * n / total:.1f}%"] for g, n in ordered]
    cell_text.append(["Total", f"{total:,}", "100.0%"])

    # bbox pins the table to fill the axes exactly (loc='center' otherwise
    # centers a fixed-height table inside a taller axes box, leaving a dead
    # gap above it under the title).
    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="left",
        colLoc="left",
        bbox=[0, 0, 1, 1],
        colWidths=[0.55, 0.22, 0.23],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.set_linewidth(0.8)
        if row == 0:
            cell.set_facecolor(HEADER_BG)
            cell.set_text_props(fontweight="bold", color=INK)
        elif row == len(cell_text):  # total row
            cell.set_facecolor(HEADER_BG)
            cell.set_text_props(fontweight="bold", color=INK)
        else:
            group_name = ordered[row - 1][0]
            cell.set_facecolor(HIGHLIGHT_BG if group_name == args.highlight else SURFACE)
            if col == 1 or col == 2:
                cell.set_text_props(ha="right")

    ax.set_title(
        f"{args.species}: Genomes per Group (n={total:,})",
        fontsize=13, fontweight="bold", pad=6, y=1.0,
    )
    footnote = f"Highlighted row = {args.highlight}, the marker-discovery target group"
    if excluded:
        footnote += f"  |  excluded: {', '.join(excluded)} ({n_excluded:,} genomes, not a real lineage)"
    fig.text(
        0.5, -0.01 / n_rows,
        footnote,
        ha="center", va="top", fontsize=8, color=INK2,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    species_slug = args.species.replace(" ", "_").replace(".", "")
    out_path = outdir / f"{species_slug}_lineage_counts_table.png"
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}", file=sys.stderr)

    # Also drop a plain TSV alongside -- easier to paste into a slide table
    # directly, or hand to someone who wants the raw numbers.
    tsv_path = outdir / f"{species_slug}_lineage_counts.tsv"
    with open(tsv_path, "w") as f:
        f.write("group\tgenomes\tpct_of_total\n")
        for g, n in ordered:
            f.write(f"{g}\t{n}\t{100 * n / total:.1f}\n")
        f.write(f"Total\t{total}\t100.0\n")
    print(f"wrote {tsv_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
