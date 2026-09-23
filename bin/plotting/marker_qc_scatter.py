#!/usr/bin/env python3
"""
Final markers -- base-pair length vs GC content, coloured by post-processing QC
verdict (pass / reject). Reads the two post-processing FASTA files directly
(headers carry length=, gc=, and, for rejects, reason=), no separate TSV needed.

QC thresholds shown as a shaded band: length >= 100 bp, GC 35-60%. A marker
passes only if it lands inside the band on both axes.

Usage: python3 marker_qc_scatter.py <pass.fasta> <reject.fasta> [-o OUTDIR] [--title TITLE]
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

BLUE = "#2a78d6"
RED = "#c0392b"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"
BAND = "#eaf2fc"

LEN_MIN = 100
GC_LO, GC_HI = 35.0, 60.0

HEADER_RE = re.compile(r"length=(\d+)\s+gc=([\d.]+)")
REASON_RE = re.compile(r"reason=([a-z_]+)")

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.size": 11,
        "font.family": "DejaVu Sans",
        "text.color": INK,
        "axes.edgecolor": INK2,
        "axes.linewidth": 0.8,
    }
)


def load_headers(path):
    """[(marker_id, length, gc, reason_or_None), ...] from a FASTA's '>' header lines."""
    rows = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = HEADER_RE.search(line)
            if not m:
                continue
            marker_id = line[1:].split()[0]
            r = REASON_RE.search(line)
            rows.append((marker_id, int(m.group(1)), float(m.group(2)), r.group(1) if r else None))
    return rows


def main():
    p = argparse.ArgumentParser(description="Final-marker length vs GC QC scatter")
    p.add_argument("pass_fasta", help="post-processing PASS fasta (length=/gc= headers)")
    p.add_argument("reject_fasta", help="post-processing REJECT fasta (length=/gc=/reason= headers)")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: current dir)")
    p.add_argument("--title", default=None, help="Figure title (default: a generic suggestion)")
    p.add_argument("--group-label", default="", help="e.g. '7PET core mode' -- appended to the default title")
    args = p.parse_args()

    passed = load_headers(args.pass_fasta)
    rejected = load_headers(args.reject_fasta)
    if not passed and not rejected:
        sys.exit("No markers found in either FASTA -- check the header format (expects 'length=N gc=NN.NN').")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 6.5))

    xmax = max([l for _, l, _, _ in passed + rejected], default=LEN_MIN) * 1.08
    ax.add_patch(Rectangle((LEN_MIN, GC_LO), xmax - LEN_MIN, GC_HI - GC_LO, facecolor=BAND, edgecolor="none", zorder=0))
    ax.axvline(LEN_MIN, color=INK2, ls="--", lw=1, zorder=1)
    ax.axhline(GC_LO, color=INK2, ls="--", lw=1, zorder=1)
    ax.axhline(GC_HI, color=INK2, ls="--", lw=1, zorder=1)

    if rejected:
        ax.scatter(
            [l for _, l, _, _ in rejected], [g for _, _, g, _ in rejected],
            s=42, facecolor=RED, edgecolor="white", linewidth=0.4, alpha=0.75, zorder=2, label=f"Reject (n={len(rejected)})",
        )
    if passed:
        ax.scatter(
            [l for _, l, _, _ in passed], [g for _, _, g, _ in passed],
            s=90, facecolor=BLUE, edgecolor="white", linewidth=0.8, zorder=3, label=f"Pass (n={len(passed)})",
        )
        for mid, l, g, _ in passed:
            ax.annotate(mid, (l, g), textcoords="offset points", xytext=(6, 4), fontsize=8, color=INK)

    ax.set_xscale("log")
    ax.set_xlabel("Marker length (bp, log scale)", fontsize=11, fontweight="bold")
    ax.set_ylabel("GC content (%)", fontsize=11, fontweight="bold")
    default_title = "Final Markers: Length vs GC Content QC" + (f" ({args.group_label})" if args.group_label else "")
    ax.set_title(args.title or default_title, fontsize=13, fontweight="bold")
    ax.grid(color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor=BAND, edgecolor="none", label=f"QC band: ≥{LEN_MIN} bp, GC {GC_LO:.0f}–{GC_HI:.0f}%"))
    labels.append(f"QC band: ≥{LEN_MIN} bp, GC {GC_LO:.0f}–{GC_HI:.0f}%")
    ax.legend(handles, labels, loc="lower right", fontsize=9, frameon=True, framealpha=0.95)

    reject_reasons = sorted({r for _, _, _, r in rejected if r})
    reason_note = ", ".join(reject_reasons) if reject_reasons else "no reason recorded"
    fig.text(
        0.5, -0.03,
        "QC: length ≥ 100 bp and GC 35–60% required for primer suitability. "
        f"{len(rejected)} of {len(rejected) + len(passed)} candidate markers rejected ({reason_note}).",
        ha="center", fontsize=8.5, color=INK2,
    )

    plt.tight_layout()
    out_path = outdir / "final_markers_qc_scatter.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
