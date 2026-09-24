#!/usr/bin/env python3

# =============================================================================
# Output file column definitions
# =============================================================================
#
# {sample_id}.gpsc_summary.tsv  — top N GPSCs by gpsc_score (N = --top, default 20)
#   GPSC                  : GPSC lineage identifier
#   avg_hit_breadth_cov   : length-weighted average breadth of coverage across
#                           all query contigs that hit this GPSC;
#                           computed as sum(bases_covered) / sum(contig_lengths)
#                           for all contigs with at least one hit to the GPSC
#   unique_matched_refs   : number of distinct reference genomes in this GPSC
#                           that received at least one hit
#   gpsc_size             : total number of reference genomes in this GPSC
#   gpsc_score            : (avg_hit_breadth_cov / gpsc_size) *
#                           sqrt(unique_matched_refs / gpsc_size)
#                           length-weighted coverage normalised by GPSC size,
#                           penalised by the sqrt fraction of GPSC refs matched
#   contig_hit_fraction   : fraction of query contigs with any hit to this GPSC
#                           (contigs_with_any_hit / total_contigs_used)
#
# {sample_id}.colours_per_contig.tsv  — one row per query contig
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

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Summarise Themisto2 pseudoalignment colour hits by GPSC lineage.")

    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
        help=(
            "Path to manifest file containing 'sample_ID' column from metadata, JSONL file output from themisto "
            "threshold-pseudoalign, Path to the query's fasta file and known GPSC label (This is an optional "
            "column to fill, metadata file provided can automatically fill this blank space)."
        ),
    )

    parser.add_argument(
        "--mapping",
        required=True,
        type=Path,
        help="Path to TSV mapping file where row index = Themisto colour ID, with columns Sample_ID and GPSC.",
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
        help=(
            "Minimum fraction of contig k-mers that must be covered for a contig to be used. "
            "Default: None (no filter applied, equivalent to 0)."
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/"),
        help=("Path to the GPSC summary output files. Default: results/."),
    )

    parser.add_argument(
        "--min-unique-refs",
        type=int,
        default=10,
        help=(
            "Minimum number of unique matched references a GPSC must have to be considered "
            "for the top prediction. Default: 10."
        ),
    )

    return parser.parse_args()


# Input functions
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
    # validate sample_ID and GPSC columns in the mapping file
    # preserve colour-row ordering, colour IDs are implicit by row. setup mapping.index
    colour2sample_mapping = pd.read_csv(mapping_path, sep="\t")

    required_cols = {"Sample_ID", "GPSC"}
    missing = required_cols - set(colour2sample_mapping.columns)

    if missing:
        raise ValueError(f"Mapping file is missing required columns: {missing}. ")

    colour2sample_mapping["Sample_ID"] = colour2sample_mapping["Sample_ID"].astype(str)
    colour2sample_mapping["GPSC"] = colour2sample_mapping["GPSC"].astype(int)

    return colour2sample_mapping


def load_metadata(metadata_path):
    sep = "," if str(metadata_path).endswith(".csv") else "\t"
    metadata = pd.read_csv(metadata_path, sep=sep, low_memory=False)

    required_cols = {"Assembly_Length", "GPSC"}
    missing = required_cols - set(metadata.columns)

    if missing:
        raise ValueError(f"Metadata file is missing required columns: {missing}.")

    return metadata


def load_fasta_lengths(fasta_path):
    # Parse contig lengths from a FASTA file.
    # If the FASTA header includes an explicit length annotation such as
    # len=1234 or length=1234, use that. Otherwise compute contig length
    # from the sequence lines that follow the header.
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    def _parse_header_length(header):
        for marker in ("len=", "length=", "LN=", "LN:"):
            if marker in header:
                start = header.index(marker) + len(marker)
                digits = []
                for ch in header[start:]:
                    if ch.isdigit():
                        digits.append(ch)
                    else:
                        break
                if digits:
                    return int("".join(digits))
        return None

    contig_lengths = {}
    current_name = None
    current_length = 0
    use_header_length = False

    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    contig_lengths[current_name] = current_length

                header = line[1:].strip()
                current_name = header.split()[0]
                header_length = _parse_header_length(header)

                if header_length is not None:
                    current_length = header_length
                    use_header_length = True
                else:
                    current_length = 0
                    use_header_length = False

                continue

            if current_name is None:
                continue

            if not use_header_length:
                current_length += len(line)

    if current_name is not None:
        contig_lengths[current_name] = current_length

    return contig_lengths


def resolve_gpsc_for_manifest(manifest, metadata):
    if "GPSC" not in manifest.columns:
        manifest = manifest.copy()
        manifest["GPSC"] = None

    gpsc_lookup = metadata.set_index("Sample_ID")["GPSC"].to_dict()

    manifest = manifest.copy()
    manifest["GPSC"] = manifest["GPSC"].astype(object)
    mask = manifest["GPSC"].isna()
    manifest.loc[mask, "GPSC"] = manifest.loc[mask, "Sample_ID"].map(gpsc_lookup)

    return manifest


def compute_contig_coverage(bases_covered_list, contig_length):
    if not bases_covered_list or not contig_length:
        return None
    return max(bases_covered_list) / contig_length


def parse_themisto_output(jsonl_path, colour2sample_mapping, contig_lengths, min_contig_coverage):
    themisto_output = Path(jsonl_path)

    if not themisto_output.exists():
        raise FileNotFoundError(f"Themisto output not found: {themisto_output}")

    gpsc_list = colour2sample_mapping["GPSC"].tolist()
    n_ref = len(gpsc_list)

    cumulative_bases_cov = defaultdict(int)  # sum of raw bases covered per GPSC across all hits
    contig_lengths_by_gpsc = defaultdict(int)  # sum of query contig lengths that hit each GPSC
    unique_colours_by_gpsc = defaultdict(set)
    contig_support = Counter()
    colours_per_contig = []
    skipped_contigs = []
    parsed_json_records = 0
    parsed_json_records_used = 0
    unmapped_colour_ids = set()

    with themisto_output.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            rec = json.loads(line)

            if "colours" not in rec:
                continue

            parsed_json_records += 1

            contig = rec.get("name", "unknown_contig")
            colours = rec.get("colours", [])
            bases_covered_list = rec.get("bases_covered", [])

            contig_length = contig_lengths.get(contig)
            coverage_fraction = compute_contig_coverage(bases_covered_list, contig_length)

            passed = (
                min_contig_coverage is None or coverage_fraction is None or coverage_fraction >= min_contig_coverage
            )

            contig_record = {
                "contig": contig,
                "n_matched_refs": len(colours),
                "best_ref_kmer_coverage": coverage_fraction,
                "passed_coverage_filter": passed,
            }

            if not passed:
                skipped_contigs.append(contig_record)
                continue

            colours_per_contig.append(contig_record)
            parsed_json_records_used += 1
            gpscs_on_this_contig = set()

            for idx, colour in enumerate(colours):
                if colour < 0 or colour >= n_ref:
                    unmapped_colour_ids.add(colour)
                    continue

                gpsc = gpsc_list[colour]
                bases_cov = bases_covered_list[idx] if idx < len(bases_covered_list) else 0
                cumulative_bases_cov[gpsc] += bases_cov
                unique_colours_by_gpsc[gpsc].add(colour)
                gpscs_on_this_contig.add(gpsc)

            for gpsc in gpscs_on_this_contig:
                contig_support[gpsc] += 1
                contig_lengths_by_gpsc[gpsc] += contig_length or 0

    return {
        "cumulative_bases_cov": cumulative_bases_cov,
        "contig_lengths_by_gpsc": contig_lengths_by_gpsc,
        "unique_colours_by_gpsc": unique_colours_by_gpsc,
        "contig_support": contig_support,
        "colours_per_contig": colours_per_contig,
        "skipped_contigs": skipped_contigs,
        "parsed_json_records": parsed_json_records,
        "parsed_json_records_used": parsed_json_records_used,
        "unmapped_colour_ids": unmapped_colour_ids,
    }


def make_summary_tables(results, colour2sample_mapping):
    gpsc_size = colour2sample_mapping["GPSC"].value_counts().to_dict()

    all_gpscs = set(gpsc_size)
    all_gpscs.update(results["cumulative_bases_cov"].keys())
    all_gpscs.update(results["unique_colours_by_gpsc"].keys())
    all_gpscs.update(results["contig_support"].keys())

    rows = []

    parsed_contigs_used = results["parsed_json_records_used"]

    for gpsc in sorted(all_gpscs, key=lambda x: str(x)):
        raw_bases = results["cumulative_bases_cov"].get(gpsc, 0)
        total_contig_len = results["contig_lengths_by_gpsc"].get(gpsc, 0)
        # length-weighted average breadth coverage across all query contigs hitting this GPSC
        avg_hit_breadth_cov = raw_bases / total_contig_len if total_contig_len else 0
        unique_hits = len(results["unique_colours_by_gpsc"].get(gpsc, set()))
        size = gpsc_size.get(gpsc, 0)
        contigs = results["contig_support"].get(gpsc, 0)

        gpsc_score = (avg_hit_breadth_cov / size) * (unique_hits / size) ** 0.5 if size else 0

        contig_hit_fraction = contigs / parsed_contigs_used if parsed_contigs_used > 0 else 0

        rows.append(
            {
                "GPSC": gpsc,
                "avg_hit_breadth_cov": avg_hit_breadth_cov,
                "unique_matched_refs": unique_hits,
                "gpsc_size": size,
                "gpsc_score": gpsc_score,
                "contig_hit_fraction": contig_hit_fraction,
            }
        )

    return pd.DataFrame(rows)


def print_top_tables(summary, top):
    cols = [
        "GPSC",
        "avg_hit_breadth_cov",
        "unique_matched_refs",
        "gpsc_size",
        "gpsc_score",
        "contig_hit_fraction",
    ]

    print("\nTop GPSCs by gpsc_score:")
    print(summary.sort_values("gpsc_score", ascending=False).head(top)[cols].to_string(index=False))


def main():
    args = parse_args()

    args.output.mkdir(parents=True, exist_ok=True)

    manifest = read_manifest(args.manifest)
    colour2sample_mapping = load_mapping(args.mapping)
    metadata = load_metadata(args.metadata)

    manifest = resolve_gpsc_for_manifest(manifest, metadata)
    query_rows = load_query_rows(manifest)

    for query in query_rows:
        sample_id = query["Sample_ID"]

        contig_lengths = load_fasta_lengths(query["fasta_path"])

        results = parse_themisto_output(
            jsonl_path=query["JSONL_path"],
            colour2sample_mapping=colour2sample_mapping,
            contig_lengths=contig_lengths,
            min_contig_coverage=args.min_contig_coverage,
        )

        summary = make_summary_tables(results, colour2sample_mapping)

        sorted_summary = summary.sort_values("gpsc_score", ascending=False)
        scores = sorted_summary["gpsc_score"]
        top_score = scores.iloc[0] if len(scores) > 0 else 0
        second_score = scores.iloc[1] if len(scores) > 1 else 0
        total_score = scores.sum()

        score_ratio = top_score / second_score if second_score > 0 else None
        top_score_fraction = top_score / total_score if total_score > 0 else 0

        eligible = sorted_summary[sorted_summary["unique_matched_refs"] >= args.min_unique_refs]
        top_gpsc = str(int(float(eligible.iloc[0]["GPSC"]))) if len(eligible) > 0 else "NA"
        match = (
            int(float(top_gpsc)) == int(float(str(query["known_gpsc"])))
            if query["known_gpsc"] is not None and top_gpsc != "NA"
            else None
        )

        contig_cols = ["contig", "n_matched_refs", "best_ref_kmer_coverage", "passed_coverage_filter"]
        colours_per_contig = pd.DataFrame(results["colours_per_contig"], columns=contig_cols)
        skipped_contigs = pd.DataFrame(results["skipped_contigs"], columns=contig_cols)

        sorted_summary.head(args.top).to_csv(args.output / f"{sample_id}.gpsc_summary.tsv", sep="\t", index=False)

        pd.DataFrame(
            [
                {
                    "Sample_ID": sample_id,
                    "predicted_GPSC": top_gpsc,
                    "known_GPSC": query["known_gpsc"],
                    "match": match,
                    "top_to_2nd_score_ratio": score_ratio,
                    "top_score_fraction": top_score_fraction,
                    "top_gpsc_score": top_score,
                }
            ]
        ).to_csv(args.output / f"{sample_id}.classification.tsv", sep="\t", index=False)

        colours_per_contig.sort_values("n_matched_refs", ascending=False).to_csv(
            args.output / f"{sample_id}.colours_per_contig.tsv", sep="\t", index=False
        )

        if len(skipped_contigs) > 0:
            skipped_contigs.to_csv(args.output / f"{sample_id}.skipped_contigs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
