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
import re
import sys
import pandas as pd
from pathlib import Path

# For labels with multiple semicolon-separated values (e.g. GPSC "1215;5"),
# --resolve-multi-value picks the smallest numeric value -- smaller is the
# current canonical label under GPSC v11 (higher numbers were merged into
# lower ones). Exception: "235;9" is a genuine biological blend of two
# distinct families (not a merge), so it's kept as a combined label "235_9"
# rather than resolved to either parent.
BLEND_EXCEPTIONS = {"235;9", "9;235"}


def resolve_multi_value_label(raw):
    if ";" not in raw:
        return raw
    parts = [p.strip() for p in raw.split(";")]
    if raw in BLEND_EXCEPTIONS:
        return "_".join(sorted(parts, key=int))
    return str(min(int(p) for p in parts))


def sanitize_id(raw):
    return re.sub(r"[^A-Za-z0-9._-]", "_", raw)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Themisto2 file-colors input from a metadata table."
    )
    parser.add_argument(
        "--metadata", required=True,
        help="Path to metadata file (CSV/TSV, optionally gzipped -- inferred from extension)."
    )
    parser.add_argument(
        "--sep", default=",",
        help="Field separator for --metadata, e.g. ',' (default) or $'\\t' for TSV."
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
    parser.add_argument(
        "--resolve-multi-value", action="store_true",
        help="For labels with multiple semicolon-separated values (e.g. GPSC "
             "'1215;5'), resolve to the smallest numeric value, keeping the "
             "known 235;9 blend as a combined '235_9' label instead of "
             "resolving to either parent. Off by default; requires numeric "
             "label values."
    )
    parser.add_argument(
        "--extra-missing-values", nargs="*", default=[],
        help="Additional literal string values (besides NaN and the built-in "
             "'unknown') to treat as a missing/unresolved label and drop, "
             "e.g. --extra-missing-values unclassifiable. Empty by default "
             "so existing behaviour is unchanged."
    )
    parser.add_argument(
        "--sanitize-sample-id", action="store_true",
        help="Before matching against assembly filenames, replace every "
             "character in --sample-col that isn't alphanumeric/'.'/'-'/'_' "
             "with '_' (matches the transform used when Pathogenwatch-style "
             "names containing '/', spaces, or parentheses were turned into "
             "filenames on disk). Off by default."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    metadata = pd.read_csv(args.metadata, sep=args.sep, low_memory=False)

    missing_cols = [c for c in [args.sample_col, args.label_col] if c not in metadata.columns]
    if missing_cols:
        sys.exit(f"Error: column(s) not found in metadata: {', '.join(missing_cols)}")

    # Treat NaN and literal sentinel strings (e.g. "unknown" for unresolved
    # sylph species calls, or "unclassifiable" for unresolved Pathogenwatch
    # typing) the same way: no usable label, so the sample is dropped rather
    # than grouped into a bogus catch-all label.
    missing_values = {"unknown", *args.extra_missing_values}
    is_missing = metadata[args.label_col].isna() | metadata[args.label_col].isin(missing_values)
    nan_label = is_missing.sum()
    metadata = metadata[~is_missing]

    n_multi_value = int(metadata[args.label_col].astype(str).str.contains(";").sum())
    if args.resolve_multi_value:
        metadata[args.label_col] = metadata[args.label_col].astype(str).apply(resolve_multi_value_label)

    # sort by label lexicographically, then by sample ID within each group
    metadata = metadata.sort_values(
        [args.label_col, args.sample_col],
        key=lambda col: col.astype(str)
    )

    sample_ids = metadata[args.sample_col].astype(str)
    if args.sanitize_sample_id:
        sample_ids = sample_ids.apply(sanitize_id)

    if args.assembly_dir:
        assembly_dir = Path(args.assembly_dir)
        assembly_files = set(os.listdir(assembly_dir))
        metadata["filename"] = sample_ids + args.assembly_suffix
        no_assembly = ~metadata["filename"].isin(assembly_files)
        on_disk_no_metadata = assembly_files - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(assembly_files)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: str(assembly_dir / f))
    else:
        with open(args.assembly_paths) as fh:
            path_list = [line.strip() for line in fh if line.strip()]
        # build lookup: basename (without suffix) -> full path
        path_lookup = {Path(p).name: p for p in path_list}
        metadata["filename"] = sample_ids + args.assembly_suffix
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
