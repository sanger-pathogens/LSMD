#!/usr/bin/env python3

"""
Design PCR primer pairs for candidate marker sequences using primer3_core.

Input: a POST_PROCESS_MARKERS filtered-markers FASTA
('{lineage_id}_filtered_markers.fasta') from post_processing_unitigs.py. Locally
GC-imbalanced windows are soft-masked (lowercased) in that FASTA. primer3 is run
with PRIMER_LOWERCASE_MASKING=1, so it will not anchor a primer's 3' end on a
masked base but a pair may still flank a masked window -- the whole marker stays
in play, per Vignesh Shetty's spec. (This replaced an earlier approach that split
each marker into clean sub-segments and designed within them, which threw away
every pair that would have bracketed a short bad window.)

All markers are batched into one Boulder-IO input and one primer3_core call
(not one invocation per marker) -- primer3_core's startup cost dwarfs the
per-record design time, so batching is the difference between seconds and
minutes for a few thousand markers.

A marker primer3 can't find a valid pair for is not a script error -- it's
written to the "no primers" summary with primer3's own PRIMER_*_EXPLAIN
diagnostics, and design continues for the rest.

Output: '{label}_primers.tsv' (one row per returned primer pair, ranked by
primer3's own pair penalty) and '{label}_no_primers.tsv'.
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

##############################################################################
# FASTA -> Boulder-IO input


@dataclass
class FastaRecord:
    id: str
    description: str  # full header line, id included, no leading '>'
    seq: str  # case preserved -- lowercase = soft-masked non-designable region


def parse_fasta(path: Path) -> Iterator[FastaRecord]:
    """Minimal FASTA reader -- no Biopython dependency, so this keeps running in
    the bare primer3_core container. Sequence case is preserved deliberately:
    the soft-mask from post_processing_unitigs.py is what primer3 keys off."""
    header: str = None
    seq_lines: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header.split()[0], header, "".join(seq_lines))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield FastaRecord(header.split()[0], header, "".join(seq_lines))


def build_boulder_input(
    fasta_path: Path,
    product_size_range: str,
    num_return: int,
    opt_size: int,
    min_size: int,
    max_size: int,
    opt_tm: float,
    min_tm: float,
    max_tm: float,
    min_gc: float,
    max_gc: float,
) -> Tuple[str, int, int]:
    """Returns (boulder_text, n_markers, n_markers_with_masking)."""
    records = []
    n_markers = 0
    n_masked = 0
    for rec in parse_fasta(fasta_path):
        n_markers += 1
        if any(c.islower() for c in rec.seq):
            n_masked += 1
        block = [
            f"SEQUENCE_ID={rec.id}",
            f"SEQUENCE_TEMPLATE={rec.seq}",
            "PRIMER_TASK=generic",
            "PRIMER_PICK_LEFT_PRIMER=1",
            "PRIMER_PICK_RIGHT_PRIMER=1",
            "PRIMER_PICK_INTERNAL_OLIGO=0",
            # soft-masked (lowercase) bases: primer3 keeps a primer's 3' end off
            # them but allows the rest of the primer, and a pair, to span them.
            "PRIMER_LOWERCASE_MASKING=1",
            f"PRIMER_NUM_RETURN={num_return}",
            f"PRIMER_PRODUCT_SIZE_RANGE={product_size_range}",
            f"PRIMER_OPT_SIZE={opt_size}",
            f"PRIMER_MIN_SIZE={min_size}",
            f"PRIMER_MAX_SIZE={max_size}",
            f"PRIMER_OPT_TM={opt_tm}",
            f"PRIMER_MIN_TM={min_tm}",
            f"PRIMER_MAX_TM={max_tm}",
            f"PRIMER_MIN_GC={min_gc}",
            f"PRIMER_MAX_GC={max_gc}",
            "=",
        ]
        records.append("\n".join(block))

    return "\n".join(records) + "\n", n_markers, n_masked


##############################################################################
# Boulder-IO output -> TSV


def parse_boulder_records(text: str) -> Iterator[Dict[str, str]]:
    """Split primer3_core's Boulder-IO stdout into one dict per '='-terminated record."""
    record: Dict[str, str] = {}
    for line in text.splitlines():
        if line == "=":
            if record:
                yield record
            record = {}
            continue
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        record[key] = value
    if record:
        yield record


def pairs_from_record(record: Dict[str, str]) -> List[Dict[str, str]]:
    """One row per PRIMER_PAIR_<i>_* pair primer3 actually returned. Positions are
    0-based and marker-relative (PRIMER_FIRST_BASE_INDEX defaults to 0), so a row
    is directly usable against the marker FASTA."""
    marker_id = record.get("SEQUENCE_ID", "")
    n_returned = int(record.get("PRIMER_PAIR_NUM_RETURNED", "0") or "0")
    rows = []
    for i in range(n_returned):
        left_pos, left_len = record[f"PRIMER_LEFT_{i}"].split(",")
        right_pos, right_len = record[f"PRIMER_RIGHT_{i}"].split(",")
        rows.append(
            {
                "marker_id": marker_id,
                "pair_rank": str(i),
                "left_seq": record[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "left_tm": record[f"PRIMER_LEFT_{i}_TM"],
                "left_pos": left_pos,
                "left_len": left_len,
                "right_seq": record[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "right_tm": record[f"PRIMER_RIGHT_{i}_TM"],
                "right_pos": right_pos,
                "right_len": right_len,
                "product_size": record[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "pair_penalty": record.get(f"PRIMER_PAIR_{i}_PENALTY", ""),
            }
        )
    return rows


TSV_COLUMNS = [
    "marker_id",
    "pair_rank",
    "left_seq",
    "left_tm",
    "left_pos",
    "left_len",
    "right_seq",
    "right_tm",
    "right_pos",
    "right_len",
    "product_size",
    "pair_penalty",
]


def write_tsv(rows: List[Dict[str, str]], path: Path):
    with open(path, "w") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(row[c] for c in TSV_COLUMNS) + "\n")


def write_no_primers(records: List[Dict[str, str]], path: Path):
    with open(path, "w") as fh:
        fh.write("marker_id\treason\n")
        for record in records:
            marker_id = record.get("SEQUENCE_ID", "")
            reason = record.get("PRIMER_PAIR_EXPLAIN", record.get("PRIMER_ERROR", "no pairs returned"))
            fh.write(f"{marker_id}\t{reason}\n")


##############################################################################
# CLI


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("fasta", type=Path, help="Filtered markers FASTA (POST_PROCESS_MARKERS output)")
    p.add_argument("--label", required=True, help="Prefix for output filenames, e.g. lineage/species id")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--primer3-bin", default="primer3_core", help="primer3_core executable (default: on PATH)")
    p.add_argument("--product-size-range", default="70-200", help="PRIMER_PRODUCT_SIZE_RANGE (default: 70-200)")
    p.add_argument("--num-return", type=int, default=3, help="PRIMER_NUM_RETURN, pairs per marker (default: 3)")
    p.add_argument("--opt-size", type=int, default=20)
    p.add_argument("--min-size", type=int, default=18)
    p.add_argument("--max-size", type=int, default=27)
    p.add_argument("--opt-tm", type=float, default=60.0)
    p.add_argument("--min-tm", type=float, default=57.0)
    p.add_argument("--max-tm", type=float, default=63.0)
    p.add_argument("--min-gc", type=float, default=20.0)
    p.add_argument("--max-gc", type=float, default=80.0)
    return p.parse_args()


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    boulder_input, n_markers, n_masked = build_boulder_input(
        args.fasta,
        product_size_range=args.product_size_range,
        num_return=args.num_return,
        opt_size=args.opt_size,
        min_size=args.min_size,
        max_size=args.max_size,
        opt_tm=args.opt_tm,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
    )
    if n_markers == 0:
        sys.exit(f"No records found in {args.fasta} -- nothing to design primers for.")

    print(
        f"Running primer3_core on {n_markers:,} marker(s) "
        f"({n_masked:,} with soft-masked regions) in {args.fasta} ...",
        file=sys.stderr,
    )
    result = subprocess.run(
        [args.primer3_bin],
        input=boulder_input,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.exit(
            f"primer3_core exited {result.returncode} -- treating as a fatal batch failure "
            f"(a single unparseable record breaks the whole Boulder-IO batch):\n{result.stderr}"
        )

    all_rows: List[Dict[str, str]] = []
    no_primer_records: List[Dict[str, str]] = []
    markers_with_primers = set()
    for record in parse_boulder_records(result.stdout):
        rows = pairs_from_record(record)
        if rows:
            all_rows.extend(rows)
            markers_with_primers.add(rows[0]["marker_id"])
        else:
            no_primer_records.append(record)

    primers_path = args.out_dir / f"{args.label}_primers.tsv"
    no_primers_path = args.out_dir / f"{args.label}_no_primers.tsv"
    write_tsv(all_rows, primers_path)
    write_no_primers(no_primer_records, no_primers_path)

    n_with_primers = len(markers_with_primers)
    print("\nResults:", file=sys.stderr)
    print(f"  Markers submitted     : {n_markers:,}", file=sys.stderr)
    print(f"  Markers with >=1 pair : {n_with_primers:,} ({100 * n_with_primers / n_markers:.1f}%)", file=sys.stderr)
    print(f"  Markers with 0 pairs  : {n_markers - n_with_primers:,}", file=sys.stderr)
    print(f"  Total primer pairs    : {len(all_rows):,}", file=sys.stderr)
    print(f"  Wrote {primers_path}", file=sys.stderr)
    print(f"  Wrote {no_primers_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
