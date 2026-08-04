#!/usr/bin/env python3

"""
Filter and rank candidate unitigs based on:
- Length >= threshold (default: 300bp)
- GC content 40-60 %
    - Sliding window GC content check - rank by closest to 50%

so gc global calculates % for one unitig
gc sliding calculates high and low GC areas of that unitig

(ranking is length first then GC balance)
Biopython
Use package SeqIO to parse FASTA format
Use package SeqUtils to calculate GC%
"""
import argparse
import sys
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO
from Bio.SeqUtils import GC

def calculate_gc(seq: str) -> float:
    """Calculate GC% for a sequence(unitig)"""
    return GC(seq)

def check_gc_window(seq:str, window_size: int, min_gc: float, max_gc: float) -> bool:
    """
    Check if GC content stays within range 40-60% across all sliding windows.

    Args:
        seq: DNA sequence
        window_size: 31 bp
        min_gc: 40%
        max_gc: 60%
    
    Returns:
        True if GC% stays in range throughout the sequence
    """
    seq = seq.upper()
    if len(seq) < window_size:
        # sequence shorter than window are dropped.
        return False

    # check every window
    for i in range(len(seq) - window_size + 1):
        window_gc = GC(seq[i:i + window_size])
        if not (min_gc <= window_gc <= max_gc):
        return False

    return True


def filter_unitigs(
    fasta_path: str,
    min_length: int = 300,
    min_gc: float = 40.0,
    max_gc: float = 60.0,
    window_size: int = 31,
) -> Tuple[List[Tuple[str, str, float, int]], List[Tuple[str, str, float, int, str]]]:
    # Filter unitigs by length and GC content.

    passed = []
    rejected = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq)
        seq_len = len(seq_str)
        gc_pct = calculate_gc(seq_str)

        # Filter by length first:
        if seq_len < min_length:
            rejected.append((record.id, seq_str, gc_pct, seq_len, f"length_too_short({seq_len}< {min_length})"))
            continue

        # Filter by global GC content second:
        if not (min_gc <= gc_pct <= max_gc):
            rejected.append((record.id, seq_str, gc_pct, seq_len, f"GC_out_of_range ({gc_pct:.2f}% not in {min_gc}-{max_gc})"))
            continue
        
        #Fo;ter by sliding window GC content (at kmer resolution) third:
        if not check_gc_window(seq_str, window_size, min_gc, max_gc):
            rejected.append((record.id, seq_str, gc_pct, seq_len, f"sliding_window_gc_failed (window size: {window_size} bp)"))
            continue

        passed.append((record.id, seq_str, gc_pct, seq_len))

    # Ranked by length - longest sequences first then by how far the GC% is from 50%
    passed.sort(key=lambda x: (-x[3], abs(x[2] - 50.0)))
        
    return passed, rejected

def write_output(results: List[Tuple[str, str, float, int]], output_path: str):
    """Write filtered unitigs to FASTA with ranking."""
    with open(output_path, "w") as f:
        for rank, (uid, seq, gc_pct, length) in enumerate(results, 1):
            f.write(f">{uid} rank={rank} length={length} gc={gc_pct:.2f}\n")
            f.write(f"{seq}\n")
 
 
def write_rejected(rejected: List[Tuple[str, str, float, int, str]], output_path: str):
    """Write rejected unitigs to FASTA with failure reason in header."""
    with open(output_path, "w") as f:
        for uid, seq, gc_pct, length, reason in rejected:
            f.write(f">{uid} length={length} gc={gc_pct:.2f} reason={reason}\n")
            f.write(f"{seq}\n")
 
 
def plot_results(results: List[Tuple[str, str, float, int]], output_path: str, 
                 min_gc: float, max_gc: float):
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
    gcs = [x[2] for x in results]
    lengths = [x[3] for x in results]
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Filtered Unitigs — {len(results)} sequences (GC: {min_gc}–{max_gc}%)", 
                 fontsize=14, fontweight="bold")
    
    # Panel 1: Length histogram
    ax = axes[0, 0]
    ax.hist(lengths, bins=20, alpha=0.7, color="steelblue", edgecolor="black")
    ax.set_xlabel("Sequence Length (bp)")
    ax.set_ylabel("Count")
    ax.set_title("Length Distribution")
    ax.grid(axis="y", alpha=0.3)
    
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
    
    # Panel 3: Ranked scatter (GC% vs length, colored by rank)
    ax = axes[1, 0]
    scatter = ax.scatter(gcs, lengths, c=ranks, cmap="viridis", s=100, alpha=0.7, 
                        edgecolors="black", linewidth=0.5)
    ax.set_xlabel("GC Content (%)")
    ax.set_ylabel("Sequence Length (bp)")
    ax.set_title("Ranked by Length → GC Distance from 50%")
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Rank")
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
    ax.set_title(f"GC%% Profile ({gc_sorted[0]:.1f}–{gc_sorted[-1]:.1f}%)")
    ax.grid(alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Saved visualization to {output_path}", file=sys.stderr)
 
 
def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter unitigs by length and GC content, with optional visualization"
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument(
        "-o", "--output", default="filtered_unitigs.fasta",
        help="Output FASTA file (default: filtered_unitigs.fasta)"
    )
    parser.add_argument(
        "-l", "--min-length", type=int, default=300,
        help="Minimum sequence length in bp (default: 300)"
    )
    parser.add_argument(
        "-g", "--gc-min", type=float, default=40.0,
        help="Minimum GC%% (default: 40.0)"
    )
    parser.add_argument(
        "-G", "--gc-max", type=float, default=60.0,
        help="Maximum GC%% (default: 60.0)"
    )
    parser.add_argument(
        "-w", "--window-size", type=int, default=31,
        help="Sliding window size for GC check in bp (default: 31, kmer size)"
    )
    parser.add_argument(
        "--plot", action="store_true",
        help="Generate visualization plots (requires matplotlib)"
    )
    parser.add_argument(
        "-p", "--plot-output", default="unitig_analysis.png",
        help="Output plot file (default: unitig_analysis.png)"
    )
    parser.add_argument(
        "-r", "--reject-output", default=None,
        help="Output FASTA file for rejected unitigs (optional)"
    )
    
    return parser.parse_args()
 
 
def main():
    args = parse_args()
    
    # Validate input file
    if not Path(args.input).exists():
        print(f"Error: {args.input} not found", file=sys.stderr)
        sys.exit(1)
    
    # Filter
    print(f"Filtering {args.input}...", file=sys.stderr)
    passed, rejected = filter_unitigs(
        args.input,
        min_length=args.min_length,
        min_gc=args.gc_min,
        max_gc=args.gc_max,
        window_size=args.window_size,
    )
    
    # Output passed
    write_output(passed, args.output)
    
    # Output rejected (if requested)
    if args.reject_output:
        write_rejected(rejected, args.reject_output)
    
    # Summary
    print(f"\nResults:", file=sys.stderr)
    print(f"  Input: {args.input}", file=sys.stderr)
    print(f"  Output (passed): {args.output}", file=sys.stderr)
    print(f"  Passed filters: {len(passed)}", file=sys.stderr)
    print(f"  Rejected: {len(rejected)}", file=sys.stderr)
    if args.reject_output:
        print(f"  Output (rejected): {args.reject_output}", file=sys.stderr)
    if passed:
        lengths = [x[3] for x in passed]
        gcs = [x[2] for x in passed]
        print(f"  Length range: {min(lengths)}–{max(lengths)} bp", file=sys.stderr)
        print(f"  GC%% range: {min(gcs):.1f}–{max(gcs):.1f}%", file=sys.stderr)
    
    # Plot (if requested)
    if args.plot:
        plot_results(passed, args.plot_output, args.gc_min, args.gc_max)
 
 
if __name__ == "__main__":
    main()
 