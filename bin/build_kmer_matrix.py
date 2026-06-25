#!/usr/bin/env python3
"""
build_kmer_matrix.py

Builds a sparse k-mer x GPSC genome-count matrix from Themisto export output,
using the color mapping produced by color_mapping_gpsc.py. Also writes one
gzipped FASTA per GPSC lineage containing all k-mers present in at least one
genome of that lineage (pre-filtering; strict filtering is done downstream by
candidate_marker_filtering.py).

Pipeline position:
    color_mapping_gpsc.py  →  [themisto build + export]  →  build_kmer_matrix.py
                                                          →  candidate_marker_filtering.py

Inputs:
    --file-colors     : file_colors_input.txt from color_mapping_gpsc.py
                        one assembly path per line; row index = Themisto color ID
    --gpsc-mapping    : gpsc_mapping.tsv from color_mapping_gpsc.py
                        columns: Sample_ID, GPSC
    --color-sets      : export.color_sets.txt from Themisto export
                        format per line: color_set_id=N size=M c1 c2 c3 ...
    --unitigs         : export.unitigs.fa from Themisto export
                        header format: > unitig_id=N color_set_id=M
    --assembly-suffix : suffix to strip from assembly filename to recover Sample_ID
    --out-dir         : output directory

Outputs:
    kmer_gpsc_matrix.tsv.gz        : sparse matrix — columns: unitig_id, gpsc, genome_count
    gpsc_<N>_kmers.fa.gz           : per-GPSC gzipped FASTA of k-mers (one per lineage)
    gpsc_genome_counts.tsv         : GPSC → total genomes in the index
    kmer_gpsc_stats.tsv            : per-GPSC k-mer count and genome count (pre-filtering)
    stats.txt                      : run summary

# NOTE for developer Jarno: consider adding a --gzip flag to `themisto export`
# so that export.unitigs.fa and export.color_sets.txt are written compressed.
# These files can be several GB; native gzip support would significantly reduce
# disk I/O overhead when running on HPC.
"""

import argparse
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build sparse k-mer x GPSC count matrix from Themisto export.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--file-colors", required=True, type=Path,
        help="file_colors_input.txt from color_mapping_gpsc.py (one assembly path per line).",
    )
    parser.add_argument(
        "--gpsc-mapping", required=True, type=Path,
        help="gpsc_mapping.tsv from color_mapping_gpsc.py (columns: Sample_ID, GPSC).",
    )
    parser.add_argument(
        "--color-sets", required=True, type=Path,
        help="export.color_sets.txt from Themisto export.",
    )
    parser.add_argument(
        "--unitigs", required=True, type=Path,
        help="export.unitigs.fa from Themisto export.",
    )
    parser.add_argument(
        "--assembly-suffix", default=".contigs.fasta",
        help="Suffix stripped from assembly filename to recover Sample_ID (default: .contigs.fasta).",
    )
    parser.add_argument(
        "--out-dir", required=True, type=Path,
        help="Output directory.",
    )
    return parser.parse_args()


def load_color_to_gpsc(file_colors_path, gpsc_mapping_path, assembly_suffix):
    """
    Build a list mapping color_id (index position) → GPSC (int).

    color_id is implicitly assigned by Themisto as the 0-based row index of
    file_colors_input.txt, so the order of paths in that file must match the
    order used when building the Themisto index.

    Returns:
        color_to_gpsc     : list[int | None], index = color_id
        gpsc_genome_counts: Counter {gpsc: n_genomes}
    """
    gpsc_df = pd.read_csv(gpsc_mapping_path, sep="\t")
    sample_to_gpsc = dict(zip(gpsc_df["Sample_ID"].astype(str), gpsc_df["GPSC"].astype(int)))

    color_to_gpsc = []
    gpsc_genome_counts = Counter()
    unmapped = 0

    with open(file_colors_path) as fh:
        for line in fh:
            path = line.strip()
            if not path:
                continue
            sample_id = Path(path).name
            if sample_id.endswith(assembly_suffix):
                sample_id = sample_id[: -len(assembly_suffix)]
            gpsc = sample_to_gpsc.get(sample_id)
            if gpsc is None:
                unmapped += 1
            else:
                gpsc_genome_counts[gpsc] += 1
            color_to_gpsc.append(gpsc)

    if unmapped:
        print(f"Warning: {unmapped} colors could not be mapped to a GPSC.", file=sys.stderr)

    return color_to_gpsc, gpsc_genome_counts


def load_unitig_data(unitigs_path):
    """
    Stream export.unitigs.fa to build:
      - colorset_to_unitigs : {color_set_id: [unitig_id, ...]}
      - unitig_to_seq       : {unitig_id: sequence}

    Both are needed: colorset_to_unitigs drives the matrix build;
    unitig_to_seq is used when writing per-GPSC FASTAs.

    Returns:
        colorset_to_unitigs: defaultdict(list)
        unitig_to_seq      : dict {unitig_id (int): sequence (str)}
        total_unitigs      : int
    """
    colorset_to_unitigs = defaultdict(list)
    unitig_to_seq = {}
    total_unitigs = 0
    current_id = None
    current_colorset = None

    with open(unitigs_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                header = line[1:].strip()
                tokens = {kv.split("=")[0]: kv.split("=")[1] for kv in header.split() if "=" in kv}
                current_id = int(tokens["unitig_id"])
                current_colorset = int(tokens["color_set_id"])
                colorset_to_unitigs[current_colorset].append(current_id)
                total_unitigs += 1
            elif current_id is not None:
                unitig_to_seq[current_id] = line

    return colorset_to_unitigs, unitig_to_seq, total_unitigs


def build_matrix(color_sets_path, color_to_gpsc, colorset_to_unitigs, out_dir):
    """
    Stream export.color_sets.txt. For each color set, compute per-GPSC genome
    counts, write sparse matrix rows, and record which unitig_ids belong to
    each GPSC (for FASTA output).

    Returns:
        kmer_gpsc_counts  : Counter {gpsc: number of k-mers present for that GPSC}
        gpsc_to_unitig_ids: defaultdict(set) {gpsc: {unitig_id, ...}}
        total_written     : number of unitigs written to the matrix
    """
    matrix_path = out_dir / "kmer_gpsc_matrix.tsv.gz"
    kmer_gpsc_counts = Counter()
    gpsc_to_unitig_ids = defaultdict(set)
    total_written = 0
    n_colors = len(color_to_gpsc)

    with gzip.open(matrix_path, "wt") as out_fh:
        out_fh.write("unitig_id\tgpsc\tgenome_count\n")

        with open(color_sets_path) as in_fh:
            for line in in_fh:
                line = line.strip()
                if not line:
                    continue

                parts = line.split()
                color_set_id = int(parts[0].split("=")[1])

                unitig_ids = colorset_to_unitigs.get(color_set_id)
                if not unitig_ids:
                    continue

                # count genomes per GPSC for this color set
                gpsc_counts = Counter()
                for token in parts[2:]:  # skip color_set_id=N and size=M tokens
                    color_id = int(token)
                    if color_id < n_colors:
                        gpsc = color_to_gpsc[color_id]
                        if gpsc is not None:
                            gpsc_counts[gpsc] += 1

                if not gpsc_counts:
                    continue

                for unitig_id in unitig_ids:
                    for gpsc, count in sorted(gpsc_counts.items()):
                        out_fh.write(f"{unitig_id}\t{gpsc}\t{count}\n")
                        kmer_gpsc_counts[gpsc] += 1
                        gpsc_to_unitig_ids[gpsc].add(unitig_id)
                    total_written += 1

    return kmer_gpsc_counts, gpsc_to_unitig_ids, total_written


def write_gpsc_fastas(gpsc_to_unitig_ids, unitig_to_seq, all_gpscs, out_dir):
    """
    Write one gzipped FASTA per GPSC lineage containing all k-mers present in
    at least one genome of that lineage (pre-filtering).
    """
    fasta_dir = out_dir / "gpsc_fastas"
    fasta_dir.mkdir(exist_ok=True)

    for gpsc in all_gpscs:
        unitig_ids = gpsc_to_unitig_ids.get(gpsc)
        if not unitig_ids:
            continue
        fasta_path = fasta_dir / f"gpsc_{gpsc}_kmers.fa.gz"
        with gzip.open(fasta_path, "wt") as fh:
            for uid in sorted(unitig_ids):
                seq = unitig_to_seq.get(uid)
                if seq:
                    fh.write(f">unitig_{uid}\n{seq}\n")


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading color → GPSC mapping...", file=sys.stderr)
    color_to_gpsc, gpsc_genome_counts = load_color_to_gpsc(
        args.file_colors, args.gpsc_mapping, args.assembly_suffix
    )
    all_gpscs = sorted(gpsc_genome_counts)
    print(f"  {len(color_to_gpsc)} colors, {len(all_gpscs)} GPSC lineages", file=sys.stderr)

    print("Loading unitig sequences and color set mapping from FASTA...", file=sys.stderr)
    colorset_to_unitigs, unitig_to_seq, total_unitigs = load_unitig_data(args.unitigs)
    print(f"  {total_unitigs:,} unitigs, {len(colorset_to_unitigs):,} unique color sets", file=sys.stderr)

    print("Streaming color sets → building matrix...", file=sys.stderr)
    kmer_gpsc_counts, gpsc_to_unitig_ids, total_written = build_matrix(
        args.color_sets, color_to_gpsc, colorset_to_unitigs, args.out_dir
    )
    print(f"  {total_written:,} unitigs written to matrix", file=sys.stderr)

    print("Writing per-GPSC k-mer FASTAs...", file=sys.stderr)
    write_gpsc_fastas(gpsc_to_unitig_ids, unitig_to_seq, all_gpscs, args.out_dir)
    print(f"  {len(all_gpscs)} FASTA files written to {args.out_dir}/gpsc_fastas/", file=sys.stderr)

    pd.DataFrame([
        {"GPSC": gpsc, "total_genomes": gpsc_genome_counts[gpsc]}
        for gpsc in all_gpscs
    ]).to_csv(args.out_dir / "gpsc_genome_counts.tsv", sep="\t", index=False)

    pd.DataFrame([
        {
            "GPSC": gpsc,
            "total_genomes": gpsc_genome_counts.get(gpsc, 0),
            "kmers_present": kmer_gpsc_counts.get(gpsc, 0),
        }
        for gpsc in all_gpscs
    ]).to_csv(args.out_dir / "kmer_gpsc_stats.tsv", sep="\t", index=False)

    lines = [
        f"Total colors (genomes) in index : {len(color_to_gpsc)}",
        f"GPSC lineages                   : {len(all_gpscs)}",
        f"Unique color sets in FASTA      : {len(colorset_to_unitigs):,}",
        f"Total unitigs                   : {total_unitigs:,}",
        f"Unitigs written to matrix       : {total_written:,}",
    ]
    with open(args.out_dir / "stats.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    for line in lines:
        print(line)


if __name__ == "__main__":
    main()
