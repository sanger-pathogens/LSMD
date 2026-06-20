#!/usr/bin/env python3

# =============================================================================
# Output file column definitions
# =============================================================================
#
# {sample_id}.gpsc_summary.tsv  — top N GPSCs by gpsc_score (N = --top, default 20)
#   GPSC                  : GPSC lineage identifier
#   sum_contig_coverage   : sum of per-contig base coverage fractions across
#                           all contigs and all matched references in this GPSC
#                           (each contig contributes bases_covered / possible_kmers,
#                           bounded at 1 per contig-reference pair)
#   unique_matched_refs   : number of distinct reference genomes in this GPSC
#                           that received at least one hit
#   gpsc_size             : total number of reference genomes in this GPSC
#   gpsc_score            : (sum_contig_coverage / gpsc_size) *
#                           sqrt(unique_matched_refs / gpsc_size)
#                           mean per-reference coverage weighted by a sqrt
#                           diversity penalty (fraction of GPSC refs matched)
#   contig_hit_fraction   : fraction of query contigs with any hit to this GPSC
#                           (contigs_with_any_hit / total_contigs_used)
#
# {sample_id}.colors_per_contig.tsv  — one row per query contig
#   contig                : contig name from the query FASTA
#   n_matched_refs        : number of reference genomes that matched this contig
#   best_ref_kmer_coverage: fraction of contig k-mers covered by the
#                           best-matching reference (max bases_covered / possible_kmers)
#   passed_coverage_filter: whether the contig passed --min-contig-coverage
#
# {sample_id}.classification.tsv  — one row per query sample
#   Sample_ID             : sample identifier from the manifest
#   predicted_GPSC        : top-ranked GPSC by gpsc_score among those with
#                           unique_matched_refs >= --min-unique-refs
#   known_GPSC            : ground-truth GPSC from the manifest (if provided)
#   match                 : True/False whether predicted matches known GPSC
#   top_to_2nd_score_ratio: gpsc_score of rank-1 / gpsc_score of rank-2;
#                           ratio close to 1.0 means the call is uncertain
#   top_score_fraction    : rank-1 gpsc_score / sum of all gpsc_scores;
#                           fraction of total evidence pointing to the top GPSC
#   top_gpsc_score        : raw gpsc_score of the predicted GPSC
# =============================================================================

import gzip
import json
import argparse
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarise Themisto2 pseudoalignment color hits by GPSC lineage."
    )

    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help="Path to manifest file containing 'sample_ID' column from metadata, JSONL file output from themisto threshold-pseudoalign, Path to the query's fasta file and known GPSC label (This is an optional column to fill, metadata file provided can automatically fill this blank space).",
    )

    parser.add_argument(
        "--mapping",
        required=True,
        type=Path,
        help="Path to TSV mapping file where row index = Themisto color ID, with columns Sample_ID and GPSC.",
    )

    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help=(
            "Path to metadata of reference genomes that make up the index the query is compared against. "
            "Required columns: 'Assembly_length' and 'GPSC'."
        ),
    )

    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Number of top GPSCs to print. Default: 20.",
    )

    parser.add_argument(
        "--min-contig-coverage",
        type=float,
        default=None,
        help="Minimum coverage filter threshold for each contig in a query. No default. Adjust to user's preference",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/"),
        help=(
            "Path to the GPSC summary output files. Default: results/."
        ),
    )

    parser.add_argument(
        "--min-unique-refs",
        type=int,
        default=10,
        help="Minimum number of unique matched references a GPSC must have to be considered for the top prediction. Default: 10.",
    )

    return parser.parse_args()

#Input functions
def read_manifest(manifest_path):
    # validate required columns in the manifest file
    # keep full dataframe for later use, including optional columns like GPSC

    manifest = pd.read_csv(manifest_path, sep="\t")
    required_cols = {"Sample_ID", "JSONL_path", "fasta_path"}

    missing = required_cols - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest file is missing required columns: {missing}.")

    return manifest


def load_query_rows(manifest):
    # Normalize and validate each row in a multi-query manifest.
    # This helper turns a raw manifest row into a predictable record,
    # so downstream code can just use sample_id, jsonl_path, fasta_path, and known_gpsc.
    query_rows = []

    for _, row in manifest.iterrows():
        query_row = {
            "Sample_ID": str(row["Sample_ID"]),
            "JSONL_path": Path(row["JSONL_path"]),
            "fasta_path": Path(row["fasta_path"]),
            # known_gpsc may be blank; keep None if not supplied
            "known_gpsc": None if pd.isna(row.get("GPSC", None)) else str(row.get("GPSC")),
        }

        # You can also normalize alternative column names here if needed.
        query_rows.append(query_row)

    return query_rows


def load_mapping(mapping_path):
    mapping = pd.read_csv(mapping_path, sep="\t")

    required_cols = {"Sample_ID", "GPSC"}
    missing = required_cols - set(mapping.columns)

    if missing:
        raise ValueError(
            f"Mapping file is missing required columns: {missing}. "
            f"Found columns: {list(mapping.columns)}"
        )

    mapping["Sample_ID"] = mapping["Sample_ID"].astype(str)
    mapping["GPSC"] = mapping["GPSC"].astype(str)

    return mapping


def parse_themisto_output(themisto_output, mapping, skip_contigs_over=None):
    raw_gpsc_hits = Counter()
    unique_colors_by_gpsc = defaultdict(set)
    contig_support = Counter()

    colors_per_contig = []
    skipped_contigs = []

    parsed_json_records = 0
    parsed_json_records_used = 0

    total_color_hits = 0
    unmapped_color_ids = set()

    with open(themisto_output) as f:
        for line in f:
            line = line.strip()

            # Skip LSF text/logs. Themisto records start with JSON.
            if not line.startswith("{"):
                continue

            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue

            if "colors" not in rec:
                continue

            parsed_json_records += 1

            contig = rec.get("name", "unknown_contig")
            colors = rec.get("colors", [])

            colors_per_contig.append(
                {
                    "contig": contig,
                    "n_colors": len(colors),
                }
            )

            if skip_contigs_over is not None and len(colors) > skip_contigs_over:
                skipped_contigs.append(
                    {
                        "contig": contig,
                        "n_colors": len(colors),
                    }
                )
                continue

            parsed_json_records_used += 1

            gpscs_on_this_contig = set()

            for color in colors:
                if color < 0 or color >= len(mapping):
                    unmapped_color_ids.add(color)
                    continue

                gpsc = mapping.iloc[color]["GPSC"]

                raw_gpsc_hits[gpsc] += 1
                unique_colors_by_gpsc[gpsc].add(color)
                gpscs_on_this_contig.add(gpsc)
                total_color_hits += 1

            for gpsc in gpscs_on_this_contig:
                contig_support[gpsc] += 1

    return {
        "raw_gpsc_hits": raw_gpsc_hits,
        "unique_colors_by_gpsc": unique_colors_by_gpsc,
        "contig_support": contig_support,
        "colors_per_contig": colors_per_contig,
        "skipped_contigs": skipped_contigs,
        "parsed_json_records": parsed_json_records,
        "parsed_json_records_used": parsed_json_records_used,
        "total_color_hits": total_color_hits,
        "unmapped_color_ids": unmapped_color_ids,
    }


def make_summary_tables(results, mapping):
    gpsc_size = mapping["GPSC"].value_counts().to_dict()

    all_gpscs = set(gpsc_size)
    all_gpscs.update(results["raw_gpsc_hits"].keys())
    all_gpscs.update(results["unique_colors_by_gpsc"].keys())
    all_gpscs.update(results["contig_support"].keys())

    rows = []

    total_color_hits = results["total_color_hits"]
    parsed_contigs_used = results["parsed_json_records_used"]

    for gpsc in sorted(all_gpscs, key=lambda x: str(x)):
        raw_hits = results["raw_gpsc_hits"].get(gpsc, 0)
        unique_hits = len(results["unique_colors_by_gpsc"].get(gpsc, set()))
        size = gpsc_size.get(gpsc, 0)
        contigs = results["contig_support"].get(gpsc, 0)

        # % of all observed color hits that belong to this GPSC.
        # Biased by GPSC size.
        raw_hit_percent = (
            100 * raw_hits / total_color_hits
            if total_color_hits
            else 0
        )

        # % of reference colors in that GPSC that were hit at least once.
        # Can become uninformative if broad contigs hit nearly everything.
        normalised_unique_fraction = (
            unique_hits / size
            if size
            else 0
        )

        # Better cross-GPSC comparison:
        # observed hits for this GPSC divided by the maximum possible hits
        # for that GPSC across the contigs that were actually used.
        normalised_raw_percent = (
            100 * raw_hits / (size * parsed_contigs_used)
            if size and parsed_contigs_used
            else 0
        )

        rows.append(
            {
                "GPSC": gpsc,
                "raw_color_hits": raw_hits,
                "raw_hit_percent": raw_hit_percent,
                "normalised_raw_percent": normalised_raw_percent,
                "unique_hit_colors": unique_hits,
                "gpsc_reference_size": size,
                "normalised_unique_fraction": normalised_unique_fraction,
                "contigs_supporting_gpsc": contigs,
            }
        )

    summary = pd.DataFrame(rows)

    return summary


def print_top_tables(summary, top):
    cols = [
        "GPSC",
        "raw_color_hits",
        "raw_hit_percent",
        "normalised_raw_percent",
        "gpsc_reference_size",
        "unique_hit_colors",
        "normalised_unique_fraction",
        "contigs_supporting_gpsc",
    ]

    print("\nTop GPSCs by NORMALISED raw percentage:")
    print(
        summary.sort_values("normalised_raw_percent", ascending=False)
        .head(top)[cols]
        .to_string(index=False)
    )

    print("\nTop GPSCs by raw hit percentage:")
    print(
        summary.sort_values("raw_hit_percent", ascending=False)
        .head(top)[cols]
        .to_string(index=False)
    )

    print("\nTop GPSCs by UNIQUE hit colors:")
    print(
        summary.sort_values("unique_hit_colors", ascending=False)
        .head(top)[cols]
        .to_string(index=False)
    )

    print("\nTop GPSCs by NORMALISED unique-color fraction:")
    print(
        summary.sort_values("normalised_unique_fraction", ascending=False)
        .head(top)[cols]
        .to_string(index=False)
    )


def main():
    args = parse_args()

    themisto_output = Path(args.themisto_output)
    mapping_path = Path(args.mapping)

    if not themisto_output.exists():
        raise FileNotFoundError(f"Themisto output not found: {themisto_output}")

    if not mapping_path.exists():
        raise FileNotFoundError(f"Mapping file not found: {mapping_path}")

    mapping = load_mapping(mapping_path)

    print("Loaded mapping file")
    print(f"Mapping path: {mapping_path}")
    print(f"Mapping rows / Themisto colors: {len(mapping)}")
    print(f"Mapping columns: {list(mapping.columns)}")

    if args.query_sample:
        matches = mapping[mapping["Sample_ID"] == args.query_sample]

        print("\nQuery sample check:")
        if len(matches) == 0:
            print(f"{args.query_sample} was NOT found exactly in the mapping.")
        else:
            for idx, row in matches.iterrows():
                print(
                    f"color_id={idx}\t"
                    f"Sample_ID={row['Sample_ID']}\t"
                    f"GPSC={row['GPSC']}"
                )

    results = parse_themisto_output(
        themisto_output=themisto_output,
        mapping=mapping,
        skip_contigs_over=args.skip_contigs_over,
    )

    print("\nParsed Themisto output")
    print(f"Themisto output: {themisto_output}")
    print(f"Parsed JSON records total: {results['parsed_json_records']}")
    print(f"Parsed JSON records used: {results['parsed_json_records_used']}")
    print(f"Total color hits counted: {results['total_color_hits']}")
    print(f"Unmapped color IDs: {len(results['unmapped_color_ids'])}")

    if args.skip_contigs_over is not None:
        print(
            f"Skipped contigs over {args.skip_contigs_over} color hits: "
            f"{len(results['skipped_contigs'])}"
        )

    colors_per_contig = pd.DataFrame(results["colors_per_contig"])
    skipped_contigs = pd.DataFrame(results["skipped_contigs"])

    summary = make_summary_tables(results, mapping)

    print_top_tables(summary, args.top)

    out_prefix = args.out_prefix

    summary_out = f"{out_prefix}.gpsc_summary.tsv"
    contig_out = f"{out_prefix}.colors_per_contig.tsv"
    skipped_out = f"{out_prefix}.skipped_contigs.tsv"

    summary.sort_values(
        [
            "normalised_raw_percent",
            "raw_hit_percent",
            "normalised_unique_fraction",
            "unique_hit_colors",
            "raw_color_hits",
        ],
        ascending=[False, False, False, False, False],
    ).to_csv(summary_out, sep="\t", index=False)

    colors_per_contig.sort_values("n_colors", ascending=False).to_csv(
        contig_out, sep="\t", index=False
    )

    if len(skipped_contigs) > 0:
        skipped_contigs.sort_values("n_colors", ascending=False).to_csv(
            skipped_out, sep="\t", index=False
        )

    print("\nSaved output files:")
    print(summary_out)
    print(contig_out)

    if len(skipped_contigs) > 0:
        print(skipped_out)

    print("\nMain interpretation:")
    print(
        "Use 'normalised_raw_percent' first for comparing GPSCs of different sizes. "
        "Use 'raw_hit_percent' only as the percentage of all observed hits. "
        "If many contigs hit tens of thousands of colors, the run is still too broad."
    )


if __name__ == "__main__":
    main()