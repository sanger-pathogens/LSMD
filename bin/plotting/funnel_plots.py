#!/usr/bin/env python3
"""
K-mer and unitig count funnel plots from CHECKPOINT's combined pipeline_counts.tsv.

Two side-by-side horizontal bar panels per target group -- k-mer count and unitig
count -- both log-x, over a fixed, curated set of three pipeline milestones (not
every checkpointed stage; this is the presentation view, not the full diagnostic
funnel):

  1. Themisto2 exported unitigs                (species_export_unitigs -- the
                                                unitigs re-exported from the built
                                                Themisto2 index; this, not the
                                                GGCAT build or Themisto2's own
                                                internal stats, is what
                                                LINEAGE_SPECIFICITY_FILTER actually
                                                consumes, so it's the true starting
                                                point of this funnel)
  2. Lineage-Specificity Filtering            (candidate_specificity_filter)
  3. Themisto2 Pseudoalignment (Markers vs ATB) (candidate_dumped_fasta's numbers --
                                                pseudoalignment scores the candidate
                                                set, it doesn't change it; the actual
                                                filtering is FILTER_ATB_MARKERS, a
                                                later, separate step)

Why estimate k-mer counts at all: only `species_themisto_index` has a REAL k-mer
count (from `themisto2 stats`), and that stage isn't even shown here -- every
stage in this funnel is plain FASTA with no k-mer count computed. Estimating
n_kmers ~= sum_bp - n_seqs*(k-1) for those (each unitig of length L contributes
L-k+1 k-mers) fills the gap, matching Jarno/Florent's guidance that k-mer count --
not unitig count -- is the right way to report marker-set size (unitig count
swings with graph recompaction even when the underlying k-mer/bp content doesn't;
see PAT-3592/PAT-3582). Estimated bars are visually distinguished (hatched, lower
alpha).

Note on species_export_unitigs vs species_ggcat_unitigs/species_themisto_index:
these three all describe "unitigs in the species-wide index" but are measured at
three different checkpoints (GGCAT's raw build, Themisto2's own internal stats,
and a re-export from the built index) and are NOT the same number -- graph
recompaction on export shifts the unitig count (4,829,653 / ~4,829,548 vs
4,996,246 exported) even though the underlying k-mer content is unchanged. This
funnel uses the export count because that's the actual input to
LINEAGE_SPECIFICITY_FILTER; it deliberately does not show a "SBWT index" bar,
since GGCAT's build and Themisto2's own stats aren't the funnel's real starting
point and including them alongside the export count invited exactly this kind of
"why don't these match" confusion.

Usage: python3 funnel_plots.py <pipeline_counts.tsv> [-k KMER_SIZE] [-o OUTDIR]
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# (real CHECKPOINT stage name, display name) -- fixed order, top to bottom.
# candidate_dumped_fasta (Themisto2 Pseudoalignment vs ATB) deliberately excluded:
# that step scores the candidate set against ATB, it doesn't change it -- the
# actual filtering on those scores happens at markers_atb_checked, below.
CURATED_STAGES = [
    ("species_export_unitigs", "Themisto2 exported unitigs"),
    ("candidate_specificity_filter", "Lineage-Specificity Filtering"),
    ("markers_atb_checked", "ATB Species-Specificity Filtering\n(PASS markers)"),
]

# Per-stage bar colour -- index stage in blue, filtering stages (lineage -> ATB)
# each get their own colour so the funnel reads as discovery (index) vs.
# filtering (candidate narrowing), not one flat block.
STAGE_COLOR = {
    "species_export_unitigs": "#2a78d6",
    "candidate_specificity_filter": "#2ca8a0",
    "markers_atb_checked": "#e0793c",
}

BLUE = "#2a78d6"
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


def load_rows(tsv_path):
    with open(tsv_path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def as_int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def unitig_count(row):
    """n_seqs (fasta kind) or n_unitigs (themisto kind) -- same concept, two tools."""
    return as_int(row.get("n_seqs")) or as_int(row.get("n_unitigs"))


def kmer_count(row, k):
    """(count, is_estimated). None if neither a real nor estimable value exists."""
    n_kmers = as_int(row.get("n_kmers"))
    if n_kmers is not None:
        return n_kmers, False
    sum_bp = as_int(row.get("sum_bp"))
    n_seqs = as_int(row.get("n_seqs"))
    if sum_bp is not None and n_seqs is not None and sum_bp > 0:
        return max(sum_bp - n_seqs * (k - 1), 0), True
    return None, False


def build_funnel(rows, group_id, species, k):
    """The four CURATED_STAGES, in that fixed order -- species-wide ones matched
    on `species`, group-specific ones (candidate_*) matched on `group_id`."""
    by_stage = {}
    for r in rows:
        if r["stage"] not in dict(CURATED_STAGES):
            continue
        owner = species if r["stage"].startswith("species_") else group_id
        if r["id"] != owner:
            continue
        by_stage[r["stage"]] = r

    funnel = []
    for stage_name, display_name in CURATED_STAGES:
        r = by_stage.get(stage_name)
        if r is None:
            continue

        # Duplication correction: each stage's own measured dupe count
        # (populated whenever CHECKPOINT could count strand pairs on an
        # actual FASTA) drives it directly -- none of this funnel's stages
        # are the species_themisto_index special case any more.
        own_dupes = as_int(r.get("n_revcomp_dupes"))
        halve_by_dupes = own_dupes not in (None, 0)

        uc = unitig_count(r)
        if uc is not None and halve_by_dupes:
            uc = round(uc / 2)

        kc, kc_est = kmer_count(r, k)
        # Halve regardless of measured vs estimated: an estimate here is
        # sum_bp/n_seqs off the SAME duplicated FASTA (both fields counted
        # every unitig twice), so it's inflated 2x exactly like the real
        # count would be -- not halving it was a bug (sum_bp/n_seqs both
        # inflated -> the estimate is inflated too, it doesn't cancel out).
        if kc is not None and halve_by_dupes:
            kc = round(kc / 2)

        funnel.append(
            {
                "stage": stage_name,
                "display": display_name,
                "unitig_count": uc,
                "unitig_halved": halve_by_dupes,
                "kmer_count": kc,
                "kmer_estimated": kc_est,
                "kmer_halved": halve_by_dupes,
            }
        )
    return funnel


def plot_panel(ax, funnel, value_key, estimated_key, halved_key, title, xlabel, show_ylabel=False):
    labels = [f["display"] for f in funnel if f[value_key] is not None]
    values = [f[value_key] for f in funnel if f[value_key] is not None]
    colors = [STAGE_COLOR.get(f["stage"], BLUE) for f in funnel if f[value_key] is not None]
    estimated = (
        [f.get(estimated_key, False) for f in funnel if f[value_key] is not None]
        if estimated_key
        else [False] * len(values)
    )

    y = list(range(len(labels)))[::-1]  # first stage at top

    for yi, val, color, est in zip(y, values, colors, estimated):
        hatch = "///" if est else None
        alpha = 0.55 if est else 0.9
        ax.barh(yi, val, color=color, alpha=alpha, hatch=hatch, edgecolor=INK2, height=0.6)
        note = " (est.)" if est else ""
        ax.text(val * 1.15, yi, f"{val:,}{note}", va="center", ha="left", fontsize=9, color=INK)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xscale("log")
    ax.set_xlabel(xlabel, fontsize=11, fontweight="bold")
    if show_ylabel:
        ax.set_ylabel("LSMD Pipeline Stage", fontsize=11, fontweight="bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def main():
    p = argparse.ArgumentParser(description="K-mer and unitig count funnel plots from pipeline_counts.tsv")
    p.add_argument("counts_tsv", help="Path to CHECKPOINT's combined pipeline_counts.tsv")
    p.add_argument("-k", "--kmer-size", type=int, default=31, help="k-mer size used to build the index (default: 31)")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: current dir)")
    args = p.parse_args()

    rows = load_rows(args.counts_tsv)
    if not rows:
        sys.exit(f"No rows in {args.counts_tsv}")

    species = rows[0]["species"]
    group_ids = sorted({r["id"] for r in rows if as_int(r["order"]) is not None and as_int(r["order"]) >= 50})
    if not group_ids:
        sys.exit("No target-group stages (order >= 50) found -- nothing to plot.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for group_id in group_ids:
        funnel = build_funnel(rows, group_id, species, args.kmer_size)
        if not funnel:
            continue

        fig, axes = plt.subplots(1, 2, figsize=(16, 0.7 * len(funnel) + 2.5))
        # TITLE TBC -- suggestion, meant to be improved:
        fig.suptitle(
            f"From Species-wide Index to ATB-Ready Candidates: {species} {group_id} Marker Discovery",
            fontsize=14,
            fontweight="bold",
        )

        plot_panel(
            axes[0],
            funnel,
            "kmer_count",
            "kmer_estimated",
            "kmer_halved",
            "K-mer count by stage",
            "K-mers Count (log scale)",
            show_ylabel=True,
        )
        plot_panel(
            axes[1], funnel, "unitig_count", None, "unitig_halved", "Unitig count by stage", "Unitigs Count (log scale)"
        )

        plt.tight_layout(rect=[0, 0.02, 1, 0.93])
        out_path = outdir / f"{species}_{group_id}_count_funnel.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
