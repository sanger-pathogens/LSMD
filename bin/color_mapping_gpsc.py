#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Themisto2 file-colors input ordered by GPSC for S. pneumoniae."
    )
    parser.add_argument(
        "--metadata", required=True,
        help="Path to metadata CSV file (must contain Sample_ID and GPSC columns)."
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
        help="Output directory for file_colors_input.txt, gpsc_mapping.tsv, and stats.txt."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    metadata = pd.read_csv(args.metadata, low_memory=False)

    for col in ["Sample_ID", "GPSC"]:
        if col not in metadata.columns:
            sys.exit(f"Error: required column '{col}' not found in metadata.")

    nan_gpsc = metadata["GPSC"].isna().sum()
    metadata = metadata.dropna(subset=["GPSC"])

    # for samples with multiple GPSC assignments (e.g. "1215;5"), pick the smallest numeric value.
    # smaller = current canonical label under v11 (higher labels were merged into lower ones).
    multi_gpsc = metadata["GPSC"].astype(str).str.contains(";")
    multi_gpsc_df = metadata.loc[multi_gpsc, ["Sample_ID", "GPSC"]].copy()
    metadata["GPSC"] = (
        metadata["GPSC"]
        .astype(str)
        .str.split(";")
        .apply(lambda parts: min(int(p.strip()) for p in parts))
    )

    metadata = metadata.sort_values(["GPSC", "Sample_ID"])

    if args.assembly_dir:
        assembly_dir = Path(args.assembly_dir)
        assembly_files = set(os.listdir(assembly_dir))
        metadata["filename"] = metadata["Sample_ID"].astype(str) + args.assembly_suffix
        no_assembly = ~metadata["filename"].isin(assembly_files)
        on_disk_no_metadata = assembly_files - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(assembly_files)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: str(assembly_dir / f))
    else:
        with open(args.assembly_paths) as fh:
            path_list = [line.strip() for line in fh if line.strip()]
        path_lookup = {Path(p).name: p for p in path_list}
        metadata["filename"] = metadata["Sample_ID"].astype(str) + args.assembly_suffix
        no_assembly = ~metadata["filename"].isin(path_lookup)
        on_disk_no_metadata = set(path_lookup.keys()) - set(metadata["filename"])
        metadata = metadata[metadata["filename"].isin(path_lookup)].reset_index(drop=True)
        metadata["file_path"] = metadata["filename"].apply(lambda f: path_lookup[f])

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata["file_path"].to_csv(
        out_dir / "file_colors_input.txt", index=False, header=False
    )

    gpsc_mapping = (
        metadata.groupby("Sample_ID", sort=False)["GPSC"]
        .apply(lambda x: ";".join(x.astype(str)))
        .reset_index()
    )
    gpsc_mapping.to_csv(out_dir / "gpsc_mapping.tsv", index=False, sep="\t")

    lines = [
        f"Samples in metadata with no assembly on disk (dropped): {no_assembly.sum()}",
        f"Assemblies on disk with no metadata entry (excluded): {len(on_disk_no_metadata)}",
        f"Total assemblies written: {len(metadata)}",
        f"Samples with GPSC labels: {len(gpsc_mapping)}",
        f"Samples with no GPSC assignment (dropped): {nan_gpsc}",
        f"Samples with multiple GPSC assignments (resolved to smallest): {len(multi_gpsc_df)}",
    ]
    if not multi_gpsc_df.empty:
        resolved = multi_gpsc_df.copy()
        resolved["GPSC_resolved"] = resolved["GPSC"].str.split(";").apply(
            lambda parts: min(int(p.strip()) for p in parts)
        )
        lines.append(resolved.to_string(index=False))

    with open(out_dir / "stats.txt", "w") as f:
        f.write("\n".join(lines) + "\n")

    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
