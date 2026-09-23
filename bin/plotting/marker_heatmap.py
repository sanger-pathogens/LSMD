#!/usr/bin/env python3
"""
Marker specificity heatmap -- shows the pipeline's actual result (which markers hit
which colour, and how strongly), not a validation report. Generic over what the
"colour" axis actually is: ATB species (the cross-species check, FILTER_ATB_MARKERS'
output) or within-species lineages (a species' own colour index, e.g. the
species-index/intersection-pseudoalign check) -- same shape of data either way, one
row per marker x colour hit. Reads two TSVs in that shared shape:
  <out>_validation.tsv       one row per marker, includes the target-group hit fraction
  <out>_group_detail.tsv     one row per marker x non-target-group hit, uncapped
(`jsonl_to_heatmap_matrix.py` produces both from a raw themisto2 pseudoalign .jsonl,
for either check.)

Design (dataviz method): sequential white->blue fill for magnitude (hit fraction),
never a diverging scale -- there's no meaningful "zero point to diverge around"
here, just "how strong is this hit." Off-target groups (anything that isn't the
target) get a red outline box around their column instead of a red fill, so the
magnitude channel (blue depth) and the status channel ("is this a problem") stay
separate and both readable, rather than conflating them into one scale.

Column selection: a full species-wide or ATB-wide colour set can run into the
thousands. This only plots the target column plus every non-target column that has
at least one hit above --min-display-frac (default 1%) for any marker in this run.
Columns below that threshold for every marker are summarised in the caption as a
count, not silently dropped without a trace.

Usage: python3 marker_heatmap.py <out>_validation.tsv <out>_group_detail.tsv \\
           <target_group> [-o OUTDIR] [--min-display-frac 1.0] [--column-label "ATB species"]
"""
import argparse
import csv
import sys
from collections import Counter
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
BLUE = "#2a78d6"
# sequential white -> blue, light mode
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", ["#fcfcfb", BLUE])

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


def load_validation(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_group_detail(path):
    if not Path(path).exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_matrix(validation_rows, detail_rows, target_group, min_display_frac, max_markers, max_columns):
    """marker_id -> {group_name: frac_pct}, plus the ordered column list
    (target group first, then non-target groups sorted by max frac desc)."""
    # A species-wide or ATB-wide run can carry hundreds of markers -- plotting
    # every row blows the figure up without bound (one real run hit 918 rows x
    # 93 cols = ~1.9 BILLION pixels, an unopenable file). Cap it, and don't cap
    # arbitrarily/by input order: keep the markers with the most off-target
    # signal -- those are the ones that actually show the specificity story
    # (a marker with zero off-target hits draws as one flat blue cell either way).
    if max_markers and len(validation_rows) > max_markers:
        off_target_hits = Counter(d["marker_id"] for d in detail_rows)
        validation_rows = sorted(
            validation_rows, key=lambda r: -off_target_hits.get(r["marker_id"], 0)
        )[:max_markers]

    matrix = {}
    marker_order = []
    for r in validation_rows:
        mid = r["marker_id"]
        marker_order.append(mid)
        tf = float(r["target_frac"]) if r.get("target_frac") not in (None, "") else 0.0
        matrix[mid] = {target_group: tf}
    kept_ids = set(marker_order)
    detail_rows = [d for d in detail_rows if d["marker_id"] in kept_ids]

    # the ATB cross-species check's own output calls this column "species"
    # (FILTER_ATB_MARKERS' *_species_detail.tsv); a within-species/lineage
    # check calls it "group" -- same shape of data, different header either way.
    group_col = "group" if detail_rows and "group" in detail_rows[0] else "species"
    col_max = {}
    for d in detail_rows:
        mid, grp, frac = d["marker_id"], d[group_col], float(d["frac"])
        if grp == "unknown":
            continue  # excluded colour (see FILTER_ATB_MARKERS summary), not a real off-target hit
        matrix.setdefault(mid, {})[grp] = frac
        col_max[grp] = max(col_max.get(grp, 0.0), frac)

    displayed = sorted((g for g, m in col_max.items() if m >= min_display_frac), key=lambda g: -col_max[g])
    hidden_count = sum(1 for g, m in col_max.items() if m < min_display_frac)
    n_columns_total = len(displayed)
    if max_columns and len(displayed) > max_columns:
        # already sorted by strongest off-target signal (col_max desc) --
        # truncating keeps the species that actually drive the story, same
        # principle as the row cap above.
        displayed = displayed[:max_columns]
        hidden_count += n_columns_total - max_columns
    columns = [target_group] + displayed
    return matrix, marker_order, columns, hidden_count, n_columns_total


def plot_heatmap(matrix, marker_order, columns, target_group, hidden_count, out_path, title, column_label,
                  row_label, border_min_frac, n_total=None, binary_fill=False):
    """Markers on x, species/lineages on y (target row pinned at top). A red
    border marks a cell only where that off-target row actually has a hit
    >= border_min_frac for that marker -- a 0%/white cell never gets one,
    target-row cells never get one either (that's the point of being on-target).

    binary_fill: skip the sequential blue gradient entirely -- target cells solid
    blue, every non-target cell plain white regardless of its actual value. Only
    honest to use when non-target values are all/almost-all exactly 0 already
    (e.g. the clean-PASS-markers view); a real gradient still matters whenever
    off-target cells carry meaningful nonzero magnitude worth distinguishing."""
    n_rows, n_cols = len(columns), len(marker_order)
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * n_cols + 2), max(4, 0.6 * n_rows + 2)))

    for yi, col in enumerate(columns):
        for xi, mid in enumerate(marker_order):
            val = matrix.get(mid, {}).get(col, 0.0)
            if binary_fill:
                facecolor = BLUE if col == target_group else "white"
            else:
                facecolor = SEQ_CMAP(val / 100.0)
            ax.add_patch(Rectangle((xi, yi), 1, 1, facecolor=facecolor, edgecolor=GRID, linewidth=0.5))
            if col != target_group and val >= border_min_frac:
                ax.add_patch(Rectangle((xi, yi), 1, 1, facecolor="none", edgecolor=RED, linewidth=1.3))

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_xticks([i + 0.5 for i in range(n_cols)])
    ax.set_xticklabels(marker_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels(columns, fontsize=8)
    for tick, col in zip(ax.get_yticklabels(), columns):
        if col == target_group:
            tick.set_fontweight("bold")
    ax.set_xlabel(row_label, fontsize=11, fontweight="bold")
    ax.set_ylabel(column_label, fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)

    if not binary_fill:
        sm = plt.cm.ScalarMappable(cmap=SEQ_CMAP, norm=plt.Normalize(0, 100))
        cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_label("K-mer hit fraction (%)", fontsize=10)

    legend_note = (
        f"Red outline = {row_label.lower()} hitting an off-target {column_label.lower()[:-1]} "
        f"at ≥{border_min_frac:g}% (not {target_group})"
    )
    if hidden_count:
        legend_note += f"  |  {hidden_count} below display threshold not shown"
    if n_total is not None and n_total > len(marker_order):
        legend_note += f"  |  showing {len(marker_order)} of {n_total} markers (most off-target hits first)"
    fig.text(0.5, -0.02, legend_note, ha="center", va="top", fontsize=8, color=INK2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Marker specificity heatmap from a validation.tsv/group_detail.tsv pair")
    p.add_argument("validation_tsv", help="<out>_validation.tsv")
    p.add_argument("group_detail_tsv", help="<out>_group_detail.tsv")
    p.add_argument("target_group", help="the on-target colour name, e.g. vibrio_cholerae or 7PET")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: current dir)")
    p.add_argument(
        "--min-display-frac",
        type=float,
        default=1.0,
        help="Only show a non-target column if >=1 marker hits it at least this %% (default 1.0)",
    )
    p.add_argument(
        "--column-label", default="ATB species", help="Y-axis label, e.g. 'ATB species' or 'Lineage' (default: 'ATB species')"
    )
    p.add_argument(
        "--row-label", default=None, help="X-axis label (default: '<target_group> Markers')"
    )
    p.add_argument(
        "--border-min-frac", type=float, default=None,
        help="Only draw the red off-target-hit border on a cell if its hit fraction is >= this %% "
        "(default: same value as --min-display-frac)",
    )
    p.add_argument("--title", default=None, help="Figure title (default: a generic suggestion)")
    p.add_argument(
        "--max-markers", type=int, default=80,
        help="Cap the number of marker rows plotted, keeping the ones with the most off-target hits "
        "(default: 80; pass 0 to disable -- not recommended past a few hundred markers, see module docstring)",
    )
    p.add_argument(
        "--max-columns", type=int, default=30,
        help="Cap the number of non-target columns plotted, keeping the strongest off-target hits "
        "(default: 30; pass 0 to disable)",
    )
    p.add_argument(
        "--verdict", choices=["all", "PASS", "FLAG", "ABSENT"], default="all",
        help="Only plot markers with this verdict (from the validation.tsv 'verdict' column); "
        "'all' (default) plots every scored marker regardless of verdict. Output filename gets "
        "a _<verdict> suffix when this isn't 'all'.",
    )
    p.add_argument(
        "--only-bordered", action="store_true",
        help="After building the matrix, drop any marker that has no off-target cell at/above "
        "--border-min-frac (i.e. no red border anywhere in its row) -- e.g. to show PASS markers "
        "that still carry some off-target signal just under the FLAG threshold, not every clean one.",
    )
    p.add_argument(
        "--out-name", default=None, help="Override the output filename (default: '<target_group>_marker_heatmap[_<verdict>].png')",
    )
    p.add_argument(
        "--first-n", type=int, default=None,
        help="Take only the first N markers in the validation.tsv's own row order (after any --verdict "
        "filter, before any off-target-signal sorting) -- e.g. a simple preview slice, not 'strongest N'.",
    )
    p.add_argument(
        "--binary-fill", action="store_true",
        help="Solid blue for the target column, plain white for every non-target cell regardless of "
        "value (no gradient, no colourbar) -- only honest when non-target values are all/nearly-all "
        "exactly 0 already; a real gradient still matters when off-target cells carry meaningful signal.",
    )
    p.add_argument(
        "--keep-ids-file", default=None,
        help="Restrict to marker_ids listed in this file (one per line) -- e.g. a pre-computed "
        "fwd/revcomp-deduplicated canonical-ID list. Applied after --verdict, before --first-n.",
    )
    args = p.parse_args()

    validation_rows = load_validation(args.validation_tsv)
    if not validation_rows:
        sys.exit(f"No rows in {args.validation_tsv}")
    if args.keep_ids_file:
        with open(args.keep_ids_file) as f:
            keep_ids = {line.strip() for line in f if line.strip()}
        validation_rows = [r for r in validation_rows if r["marker_id"] in keep_ids]
        if not validation_rows:
            sys.exit(f"None of the IDs in {args.keep_ids_file} found in {args.validation_tsv}.")
    if args.verdict != "all":
        validation_rows = [r for r in validation_rows if r.get("verdict") == args.verdict]
        if not validation_rows:
            sys.exit(f"No markers with verdict={args.verdict} in {args.validation_tsv} -- nothing to plot.")
    if args.first_n:
        validation_rows = validation_rows[: args.first_n]
    detail_rows = load_group_detail(args.group_detail_tsv)

    HARD_ROW_CEILING = 1000  # past this, even --max-markers 0 refuses -- protects against an unopenable multi-GB file
    if not args.max_markers and len(validation_rows) > HARD_ROW_CEILING:
        sys.exit(
            f"{len(validation_rows)} markers with --max-markers 0 -- refusing (past {HARD_ROW_CEILING} rows "
            f"the figure becomes an unopenable multi-billion-pixel file). Set an explicit --max-markers instead."
        )

    matrix, marker_order, columns, hidden_count, n_columns_total = build_matrix(
        validation_rows, detail_rows, args.target_group, args.min_display_frac, args.max_markers, args.max_columns
    )
    n_before_border_filter = len(marker_order)
    border_min_frac = args.border_min_frac if args.border_min_frac is not None else args.min_display_frac

    if args.only_bordered:
        marker_order = [
            mid for mid in marker_order
            if any(v >= border_min_frac for col, v in matrix.get(mid, {}).items() if col != args.target_group)
        ]
        if not marker_order:
            sys.exit(f"No markers have an off-target hit >= {border_min_frac:g}% -- nothing to plot with --only-bordered.")

    n_shown = len(marker_order)
    n_total = len(validation_rows) if args.max_markers else n_before_border_filter

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.verdict}" if args.verdict != "all" else ""
    if args.only_bordered:
        suffix += "_bordered"
    out_path = outdir / (args.out_name or f"{args.target_group}_marker_heatmap{suffix}.png")
    # TITLE TBC -- suggestion, meant to be improved:
    default_title = f"Candidate Markers vs {args.column_label}: {args.target_group} Specificity"
    if args.verdict != "all":
        default_title += f" ({args.verdict} only)"
    title = args.title or default_title
    row_label = args.row_label or f"{args.target_group} Markers"
    plot_heatmap(
        matrix, marker_order, columns, args.target_group, hidden_count, out_path, title,
        args.column_label, row_label, border_min_frac, n_total, binary_fill=args.binary_fill,
    )
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
