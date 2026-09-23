#!/usr/bin/env python3

"""
Filter and rank candidate unitigs based on:
- Length >= threshold (default: 100bp)
- Global GC content 35-60 %
- Sliding-window GC (31bp, kmer resolution) is not a reject -- out-of-range
  windows are soft-masked (lowercased) in the output FASTA instead, and their
  coordinates listed in the header. The rest of the unitig stays uppercase and
  usable; primer3 (run with PRIMER_LOWERCASE_MASKING=1) then keeps primer 3'
  ends off the masked bases while still allowing a pair to flank them.
  Per Vignesh Shetty.

(ranking is length first then GC balance)
Biopython
Use package SeqIO to parse FASTA format
Use package SeqUtils to calculate GC%
"""

import argparse
import statistics
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


@dataclass
class UnitigResult:
    """One unitig's filtering result -- shared by both passed and rejected lists."""

    id: str
    seq: str
    gc_pct: float
    length: int
    non_designable: List[Tuple[int, int]] = field(default_factory=list)  # passed only
    reason: Optional[str] = None  # rejected only


def calculate_gc(seq: str) -> float:
    """Calculate GC% for a sequence(unitig)"""
    return gc_fraction(seq) * 100


def find_non_designable_windows(seq: str, window_size: int, min_gc: float, max_gc: float) -> List[Tuple[int, int]]:
    """
    Find sliding-window coordinates whose GC% falls outside range -- these get
    soft-masked, not used to reject the whole unitig.

    Precondition: len(seq) >= window_size (caller's job to filter length first).

    Returns: merged list of (start, end) 0-based half-open coordinate pairs.
    """
    seq = seq.upper()
    if len(seq) < window_size:
        raise ValueError(
            f"find_non_designable_windows precondition violated: sequence length "
            f"({len(seq)}) is shorter than window_size ({window_size})."
        )

    regions = []
    for i in range(len(seq) - window_size + 1):
        window_end = i + window_size
        window_gc = gc_fraction(seq[i:window_end]) * 100
        if not (min_gc <= window_gc <= max_gc):
            regions.append((i, window_end))

    return merge_windows(regions)


def merge_windows(regions: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Collapse the overlapping 1bp-stepped windows into contiguous runs.

    find_non_designable_windows appends one (start, end) per failing position, so
    a single bad stretch comes out as ~window_size overlapping pairs. Merge them
    so the header and the soft-mask both work off a handful of real regions.
    Input is already sorted by start (scan order).
    """
    merged: List[Tuple[int, int]] = []
    for start, end in regions:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def soft_mask(seq: str, regions: List[Tuple[int, int]]) -> str:
    """Lowercase every base inside a non-designable region; the rest stays upper.

    primer3 with PRIMER_LOWERCASE_MASKING=1 won't anchor a primer's 3' end on a
    lowercase base but will still place a pair that flanks the region, so the
    whole fragment stays usable -- per Vignesh Shetty's spec.
    """
    chars = list(seq.upper())
    n = len(chars)
    for start, end in regions:
        for i in range(start, min(end, n)):
            chars[i] = chars[i].lower()
    return "".join(chars)


def dedupe_strand_pairs(records: List) -> Tuple[List, Dict[str, List[str]]]:
    """Collapse forward/reverse-complement duplicate records to one per locus.

    A de Bruijn graph is inherently double-stranded, and some dump tools (e.g.
    `sbwt dump-unitigs` -- see PAT-3592) report both strand-walks of the same
    unitig as separate records, not independent markers. Keep one
    representative per canonical sequence (whichever orientation sorts
    first); PCR primer design is strand-symmetric anyway (primer3 already
    searches both strands of whatever template it's given), so nothing is
    lost by dropping the duplicate -- only redundant downstream computation
    and doubled marker counts.

    Returns (deduped_records, dup_of) where dup_of maps the kept record's id
    to the ids of every duplicate collapsed into it (a kept record can have
    more than one duplicate, e.g. 3+ strand-walks of the same locus).
    """
    canonical_seen: Dict[str, str] = {}  # canonical seq -> kept record's id
    dup_of: Dict[str, List[str]] = defaultdict(list)  # kept record's id -> collapsed duplicates' ids
    deduped = []
    for record in records:
        seq_str = str(record.seq)
        canonical = min(seq_str, str(Seq(seq_str).reverse_complement()))
        if canonical in canonical_seen:
            dup_of[canonical_seen[canonical]].append(record.id)
            continue
        canonical_seen[canonical] = record.id
        deduped.append(record)
    return deduped, dup_of


def filter_unitigs(
    fasta_path: str,
    min_length: int = 100,
    min_gc: float = 35.0,
    max_gc: float = 60.0,
    window_size: int = 31,
) -> Tuple[List[UnitigResult], List[UnitigResult], int]:
    # Filter unitigs by length and global GC content; soft-mask (don't reject)
    # local sliding-window GC dips as non-designable regions.

    # find_non_designable_windows requires len(seq) >= window_size -- guaranteed
    # by the length filter below only if window_size <= min_length.
    if window_size > min_length:
        raise ValueError(f"window_size ({window_size}) cannot exceed min_length ({min_length})")

    records, dup_of = dedupe_strand_pairs(list(SeqIO.parse(fasta_path, "fasta")))

    passed = []
    rejected = []

    for record in records:
        seq_str = str(record.seq)
        seq_len = len(seq_str)
        gc_pct = calculate_gc(seq_str)

        # Filter by length first:
        if seq_len < min_length:
            rejected.append(
                UnitigResult(record.id, seq_str, gc_pct, seq_len, reason=f"length_too_short ({seq_len} < {min_length})")
            )
            continue

        # Filter by global GC content second:
        if not (min_gc <= gc_pct <= max_gc):
            rejected.append(
                UnitigResult(
                    record.id,
                    seq_str,
                    gc_pct,
                    seq_len,
                    reason=f"GC_out_of_range ({gc_pct:.2f}% not in {min_gc}-{max_gc})",
                )
            )
            continue

        # Sliding-window GC third: record non-designable regions (soft-masked
        # in write_output), don't reject.
        non_designable = find_non_designable_windows(seq_str, window_size, min_gc, max_gc)

        passed.append(UnitigResult(record.id, seq_str, gc_pct, seq_len, non_designable=non_designable))

    # Ranked by length - longest sequences first then by how far the GC% is from 50%
    passed.sort(key=lambda r: (-r.length, abs(r.gc_pct - 50.0)))

    return passed, rejected, sum(len(dupes) for dupes in dup_of.values())


def write_output(results: List[UnitigResult], output_path: str):
    """Write filtered unitigs to FASTA: sequence soft-masked (non-designable
    windows lowercased), ranking + non-designable coordinates in the header."""
    with open(output_path, "w") as f:
        for rank, r in enumerate(results, 1):
            regions = ",".join(f"{s}-{e}" for s, e in r.non_designable) or "none"
            masked_bp = sum(e - s for s, e in r.non_designable)
            f.write(
                f">{r.id} rank={rank} length={r.length} gc={r.gc_pct:.2f} "
                f"masked_bp={masked_bp} non_designable={regions}\n"
            )
            f.write(f"{soft_mask(r.seq, r.non_designable)}\n")


def write_rejected(rejected: List[UnitigResult], output_path: str):
    """Write rejected unitigs to FASTA with failure reason in header."""
    with open(output_path, "w") as f:
        for r in rejected:
            f.write(f">{r.id} length={r.length} gc={r.gc_pct:.2f} reason={r.reason}\n")
            f.write(f"{r.seq}\n")


def plot_results(results: List[UnitigResult], output_path: str, min_gc: float, max_gc: float):
    """
    Generate 4-panel visualization:
    - Length histogram
    - GC% histogram
    - Ranked scatter (GC% vs length)
    - GC% line plot (sorted, like Biopython tutorial)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: matplotlib required for plots. Install: pip install matplotlib", file=sys.stderr)
        return

    if not results:
        print("No data to plot", file=sys.stderr)
        return

    ranks = list(range(1, len(results) + 1))
    gcs = [r.gc_pct for r in results]
    lengths = [r.length for r in results]

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Filtered Unitigs — {len(results)} sequences (GC: {min_gc}–{max_gc}%)", fontsize=14, fontweight="bold"
    )

    # Panel 1: Length histogram
    ax = axes[0, 0]
    ax.hist(lengths, bins=20, alpha=0.7, color="steelblue", edgecolor="black")
    median_len = statistics.median(lengths)
    ax.axvline(median_len, color="orange", linestyle="--", linewidth=2, label=f"Median: {median_len:.0f} bp")
    ax.set_xlabel("Sequence Length (bp)")
    ax.set_ylabel("Count")
    ax.set_title(f"Length Distribution ({min(lengths)}–{max(lengths)} bp)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    # Panel 2: GC% histogram
    ax = axes[0, 1]
    ax.hist(gcs, bins=20, alpha=0.7, color="steelblue", edgecolor="black")
    ax.axvspan(min_gc, max_gc, alpha=0.2, color="green", label=f"Target: {min_gc}–{max_gc}%")
    ax.axvline(50.0, color="orange", linestyle="--", linewidth=2, label="Optimal (50%)")
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Count")
    ax.set_title("GC Distribution")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    # Panel 3: length vs GC% -- this is the "final answer" panel (which markers
    # survived, and are they well-centered in spec), so it gets the same target-
    # range shading Panel 2 has (was missing here before -- inconsistent) and,
    # for a small enough marker set to actually read, a direct rank label on
    # every point instead of a colorbar. A colorbar-by-rank forces a
    # cross-reference for each point; a direct label doesn't, and at the scale
    # a final marker panel is usually at (single digits to a few dozen) there's
    # room to show it. Falls back to the colorbar beyond that so it doesn't
    # become a pile of overlapping numbers on a big candidate set.
    ax = axes[1, 0]
    ax.axvspan(min_gc, max_gc, alpha=0.15, color="green", zorder=0)
    DIRECT_LABEL_MAX = 40
    if len(results) <= DIRECT_LABEL_MAX:
        ax.scatter(gcs, lengths, s=110, color="#2a78d6", alpha=0.85, edgecolors="black", linewidth=0.6, zorder=2)
        for rank, gc, length in zip(ranks, gcs, lengths):
            ax.annotate(str(rank), (gc, length), xytext=(4, 4), textcoords="offset points", fontsize=8, color="#0b0b0b")
    else:
        scatter = ax.scatter(
            gcs, lengths, c=ranks, cmap="viridis", s=80, alpha=0.7, edgecolors="black", linewidth=0.5, zorder=2
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Rank")
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Sequence Length (bp)")
    ax.set_title(f"Final markers: length vs GC% (target {min_gc}–{max_gc}%)")
    ax.axvline(50.0, color="orange", linestyle="--", linewidth=1, alpha=0.5)
    ax.grid(alpha=0.3)

    # Panel 4: GC% line plot (sorted, Biopython tutorial style)
    ax = axes[1, 1]
    gc_sorted = sorted(gcs)
    ax.plot(gc_sorted, marker="o", markersize=4, alpha=0.7, color="steelblue", linewidth=1)
    ax.axhline(50.0, color="orange", linestyle="--", linewidth=2, label="Optimal (50%)")
    ax.axhspan(min_gc, max_gc, alpha=0.1, color="green")
    ax.set_xlabel("Unitig Index (sorted by GC%)")
    ax.set_ylabel("GC Content (%)")
    ax.set_title(f"GC% Profile ({gc_sorted[0]:.1f}–{gc_sorted[-1]:.1f}%)")
    ax.grid(alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved visualization to {output_path}", file=sys.stderr)


OUTPUT_FILENAME = "filtered_unitigs.fasta"
REJECTED_FILENAME = "rejected_unitigs.fasta"
PLOT_FILENAME = "unitig_analysis.png"


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Filter unitigs by length and GC content, with optional visualization")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        help=f"Output directory (default: '.'). Filenames are fixed: {OUTPUT_FILENAME}, {REJECTED_FILENAME}, "
        f"{PLOT_FILENAME}",
    )
    parser.add_argument(
        "-l", "--min-length", type=int, default=100, help="Minimum sequence length in bp (default: 100)"
    )
    parser.add_argument("-g", "--gc-min", type=float, default=35.0, help="Minimum GC%% (default: 35.0)")
    parser.add_argument("-G", "--gc-max", type=float, default=60.0, help="Maximum GC%% (default: 60.0)")
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        default=31,
        help="Sliding window size for GC check in bp (default: 31, kmer size)",
    )
    parser.add_argument("--plot", action="store_true", help="Generate visualization plots (requires matplotlib)")
    parser.add_argument("-r", "--write-rejected", action="store_true", help="Also write rejected unitigs to FASTA")

    return parser.parse_args()


def main():
    args = parse_args()

    # Validate input file
    if not Path(args.input).exists():
        print(f"Error: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / OUTPUT_FILENAME
    reject_output_path = outdir / REJECTED_FILENAME
    plot_output_path = outdir / PLOT_FILENAME

    # Filter
    print(f"Filtering {args.input}...", file=sys.stderr)
    passed, rejected, n_dupes = filter_unitigs(
        args.input,
        min_length=args.min_length,
        min_gc=args.gc_min,
        max_gc=args.gc_max,
        window_size=args.window_size,
    )

    # Output passed
    write_output(passed, output_path)

    # Output rejected (if requested)
    if args.write_rejected:
        write_rejected(rejected, reject_output_path)

    # Summary
    print("\nResults:", file=sys.stderr)
    print(f"  Input: {args.input}", file=sys.stderr)
    print(f"  Output (passed): {output_path}", file=sys.stderr)
    if args.write_rejected:
        print(f"  Output (rejected): {reject_output_path}", file=sys.stderr)
    if n_dupes:
        print(f"  Forward/reverse-complement duplicates collapsed: {n_dupes}", file=sys.stderr)
    print(f"  Passed filters: {len(passed)}", file=sys.stderr)
    print(f"  Rejected: {len(rejected)}", file=sys.stderr)
    if passed:
        lengths = [r.length for r in passed]
        gcs = [r.gc_pct for r in passed]
        print(
            f"  Length range: {min(lengths)}–{max(lengths)} bp (median {statistics.median(lengths):.0f} bp)",
            file=sys.stderr,
        )
        print(f"  GC% range: {min(gcs):.1f}–{max(gcs):.1f}%", file=sys.stderr)

    # Plot (if requested)
    if args.plot:
        plot_results(passed, plot_output_path, args.gc_min, args.gc_max)


if __name__ == "__main__":
    main()
