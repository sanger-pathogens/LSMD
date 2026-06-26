#!/usr/bin/env python3
"""
Filters candidate marker k-mers from a Themisto export by lineage specificity.

Works for any species and grouping column (e.g. GPSC, strain, serotype) via
--label-col. Loads the sparse k-mer x group matrix once, applies two strict
mode filters across all groups in a single run, then streams the unitigs FASTA
once to collect sequences for all passing k-mers before writing per-group output.

Strict mode filters:
    Step 1: k-mer present in ALL genomes of the target group.
    Step 2: k-mer absent from ALL other groups.
    Step 3 (TODO): k-mer absent from non-target species in a GTDB LexiMap index.
"""

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter candidate marker k-mers by lineage specificity (strict mode).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--matrix", required=True, type=Path,
        help="Sparse k-mer matrix TSV.gz from build_kmer_matrix.py.",
    )
    parser.add_argument(
        "--genome-counts", required=True, type=Path,
        help="Genome counts TSV from build_kmer_matrix.py (columns: <label-col>, total_genomes).",
    )
    parser.add_argument(
        "--unitigs", required=True, type=Path,
        help="export.unitigs.fa from Themisto export; streamed once to collect sequences.",
    )
    parser.add_argument(
        "--label-col", default="gpsc",
        help="Column name for the lineage grouping in the matrix and genome-counts files (default: gpsc).",
    )
    parser.add_argument(
        "--out-dir", required=True, type=Path,
        help="Output directory for strict marker FASTAs and stats.",
    )
    # TODO: uncomment once LexiMap GTDB index is available
    # parser.add_argument(
    #     "--lexicmap-output", type=Path,
    #     help="LexicMap output against GTDB index for cross-species specificity check (step 3).",
    # )
    return parser.parse_args()


def compute_strict_ids(matrix_path, label_col, total_genomes):
    """
    Load the sparse matrix once and compute strict marker unitig_ids per group.

    Returns:
        group_strict_ids : dict {group: set of unitig_ids passing steps 1 and 2}
        group_step1_ids  : dict {group: set of unitig_ids passing step 1 only}
    """
    print("Loading matrix...", file=sys.stderr)
    matrix = pd.read_csv(matrix_path, sep="\t", compression="gzip")
    print(f"  {len(matrix):,} rows loaded", file=sys.stderr)

    if label_col not in matrix.columns:
        sys.exit(f"Error: column '{label_col}' not found in matrix. Columns: {list(matrix.columns)}")

    # step 2: unitigs appearing in exactly one group
    group_counts_per_unitig = matrix.groupby("unitig_id")[label_col].nunique()
    single_group_unitigs = set(group_counts_per_unitig[group_counts_per_unitig == 1].index)
    print(f"  {len(single_group_unitigs):,} unitigs appear in exactly one group", file=sys.stderr)

    group_strict_ids = {}
    group_step1_ids = {}

    for group, grp in matrix.groupby(label_col):
        total = total_genomes.get(group)
        if total is None:
            continue
        # step 1: present in ALL genomes of this group
        step1 = set(grp[grp["genome_count"] == total]["unitig_id"])
        # step 2: absent from all other groups
        strict = step1 & single_group_unitigs
        group_step1_ids[group] = step1
        group_strict_ids[group] = strict

    return group_strict_ids, group_step1_ids


def collect_sequences(unitigs_path, needed_ids):
    """
    Stream export.unitigs.fa once and collect sequences only for unitig_ids
    in needed_ids (the union of all strict marker IDs across all groups).

    Returns:
        sequences: dict {unitig_id (int): sequence (str)}
    """
    print("Collecting sequences from unitigs FASTA...", file=sys.stderr)
    sequences = {}
    current_id = None

    with open(unitigs_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                header = line[1:].strip()
                tokens = {kv.split("=")[0]: kv.split("=")[1] for kv in header.split() if "=" in kv}
                current_id = int(tokens["unitig_id"])
            elif current_id is not None and current_id in needed_ids:
                sequences[current_id] = line

    print(f"  {len(sequences):,} sequences collected", file=sys.stderr)
    return sequences


def write_outputs(group_strict_ids, group_step1_ids, sequences, total_genomes, all_groups, label_col, out_dir):
    """
    Write per-group gzipped strict marker FASTAs and the summary stats TSV.
    """
    all_stats = []

    for group in all_groups:
        strict_ids = group_strict_ids.get(group, set())
        step1_ids = group_step1_ids.get(group, set())
        total = total_genomes[group]

        out_fasta = out_dir / f"{label_col}_{group}_strict_markers.fa.gz"
        n_written = 0
        with gzip.open(out_fasta, "wt") as fh:
            for uid in sorted(strict_ids):
                seq = sequences.get(uid)
                if seq:
                    fh.write(f">unitig_{uid}\n{seq}\n")
                    n_written += 1

        all_stats.append({
            label_col: group,
            "total_genomes": total,
            "step1_pass_all_genomes": len(step1_ids),
            "step1_dropped": len(step1_ids) - len(strict_ids),
            "step2_pass_single_group": n_written,
            "step2_dropped": len(strict_ids) - n_written,
        })

        print(
            f"  {label_col} {group}: {len(step1_ids):,} step1 -> {n_written:,} strict markers",
            file=sys.stderr,
        )

    pd.DataFrame(all_stats).to_csv(out_dir / "strict_marker_stats.tsv", sep="\t", index=False)


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    gc_df = pd.read_csv(args.genome_counts, sep="\t")

    if args.label_col not in gc_df.columns:
        sys.exit(f"Error: column '{args.label_col}' not found in genome-counts file. Columns: {list(gc_df.columns)}")

    total_genomes = dict(zip(gc_df[args.label_col], gc_df["total_genomes"]))
    all_groups = sorted(total_genomes)
    print(f"{len(all_groups)} groups to process (label: {args.label_col})", file=sys.stderr)

    group_strict_ids, group_step1_ids = compute_strict_ids(args.matrix, args.label_col, total_genomes)

    needed_ids = set().union(*group_strict_ids.values()) if group_strict_ids else set()
    sequences = collect_sequences(args.unitigs, needed_ids)

    print("Writing per-group strict marker FASTAs...", file=sys.stderr)
    write_outputs(
        group_strict_ids, group_step1_ids, sequences,
        total_genomes, all_groups, args.label_col, args.out_dir,
    )

    print(f"\nDone. Stats written to {args.out_dir}/strict_marker_stats.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
