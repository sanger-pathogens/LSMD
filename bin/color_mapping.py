#!/usr/bin/env python3
"""
Generic color mapping script for Themisto2 --file-colors input.

Reads a metadata CSV, sorts samples by a user-specified label column
(lexicographically), matches them to assembly files on disk, and writes:
  - file_colors_input.txt  : ordered assembly paths for --file-colors
  - label_mapping.tsv      : Sample_ID -> label
  - stats.txt              : reconciliation summary
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Themisto2 file-colors input from a metadata table."
    )
    parser.add_argument(
        "--metadata", required=True,
        help="Path to metadata CSV file."
    )
    parser.add_argument(
        "--sample-col", default="Sample_ID",
        help="Column name in metadata containing sample identifiers (default: Sample_ID)."
    )
    parser.add_argument(
        "--label-col", required=True,
        help="Column name to sort samples by (lexicographically) and write to label_mapping.tsv."
    )
    asm = parser.add_mutually_exclusive_group(required=True)
    asm.add_argument(
        "--assembly-dir",
        help="Directory containing assembly FASTA files."
    )
    asm.add_argument(
        "--assembly-paths",
        help="Text file listing one assembly path per line."
    )
    parser.add_argument(
        "--assembly-suffix", default=".contigs.fasta",
        help="Suffix appended to Sample_ID to form the assembly filename (default: .contigs.fasta). Only used with --assembly-dir."
    )
    parser.add_argument(
        "--out-dir", required=True,
        help="Output directory for file_colors_input.txt, label_mapping.tsv, and stats.txt."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    metadata = pd.read_csv(args.metadata, low_memory=False)

    missing_cols = [c for c in [args.sample_col, args.label_col] if c not in metadata.columns]
    if missing_cols:
        sys.exit(f"Error: column(s) not found in metadata: {', '.join(missing_cols)}")

    nan_label = metadata[args.label_col].isna().sum()
    metadata = metadata.dropna(subset=[args.label_col])

    # sort by label lexicographically, then by sample ID within each group
    metadata = metadata.sort_values(
        [args.label_col, args.sample_col],
        key=lambda col: col.astype(str)
    )

    if args.assembly_dir:
        assembly_dir = Path(args.assembly_dir)
        assembly_files = set(os.listdir(assembly_dir))
        metadata["filename"] = metadata[args.sample_col].astype(str) + args.assembly_suffix
        no_assembly = ~metadata["filename"].isin(assembly_files)
        on_disk_no_metadata = assembly_files - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(assembly_files)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: str(assembly_dir / f))
    else:
        with open(args.assembly_paths) as fh:
            path_list = [line.strip() for line in fh if line.strip()]
        # build lookup: basename (without suffix) -> full path
        path_lookup = {Path(p).name: p for p in path_list}
        metadata["filename"] = metadata[args.sample_col].astype(str) + args.assembly_suffix
        no_assembly = ~metadata["filename"].isin(path_lookup)
        listed_names = set(path_lookup.keys())
        on_disk_no_metadata = listed_names - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(path_lookup)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: path_lookup[f])

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata["file_path"].to_csv(
        out_dir / "file_colors_input.txt", index=False, header=False
    )

    label_mapping = (
        metadata[[args.sample_col, args.label_col]]
        .drop_duplicates()
        .rename(columns={args.sample_col: "Sample_ID", args.label_col: "label"})
    )
    label_mapping.to_csv(out_dir / "label_mapping.tsv", index=False, sep="\t")

    lines = [
        f"Metadata column used for ordering: {args.label_col}",
        f"Samples with no label (dropped): {nan_label}",
        f"Samples in metadata with no assembly on disk (dropped): {no_assembly.sum()}",
        f"Assemblies on disk with no metadata entry (excluded): {len(on_disk_no_metadata)}",
        f"Total assemblies written: {len(metadata)}",
    ]

    stats_path = out_dir / "stats.txt"
    with open(stats_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
