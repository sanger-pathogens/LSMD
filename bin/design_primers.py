#!/usr/bin/env python3

"""
Design PCR primer pairs for candidate marker sequences using primer3_core.

Input: a POST_PROCESS_MARKERS filtered-markers FASTA
('{lineage_id}_filtered_markers.fasta') -- headers carry
'non_designable=<start>-<end>,...' coordinates from post_processing_unitigs.py's
sliding-window GC check. Those regions are passed to primer3 as
SEQUENCE_EXCLUDED_REGION so no primer is ever placed across a locally
GC-imbalanced stretch, without re-deriving anything post_processing_unitigs.py
already worked out.

All markers are batched into one Boulder-IO input and one primer3_core call
(not one process invocation per marker) -- primer3_core's own startup cost
dwarfs the per-record design time, so batching is the difference between
seconds and minutes for a few thousand markers.

A marker primer3 can't find a valid pair for (SEQUENCE_ID present in the
output with no PRIMER_LEFT_0_SEQUENCE, or a nonzero PRIMER_ERROR field) is
not a script error -- it's written to the "no primers" summary, and design
continues for the rest.

Output: '{label}_primers.tsv' (one row per returned primer pair, ranked by
primer3's own pair penalty) and '{label}_no_primers.tsv' (markers primer3
returned zero pairs for, with primer3's own PRIMER_*_EXPLAIN diagnostics).
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Tuple

##############################################################################
### FASTA -> Boulder-IO input


@dataclass
class FastaRecord:
    id: str
    description: str  # full header line, id included, no leading '>'
    seq: str


def parse_fasta(path: Path) -> Iterator[FastaRecord]:
    """Minimal FASTA reader -- the only thing needed here is id/header/sequence,
    so no Biopython dependency (keeps this script running in the bare
    primer3_core container, no second package/install needed at runtime)."""
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


def parse_non_designable(field: str) -> List[Tuple[int, int]]:
    """'123-154,301-332' -> [(123, 154), (301, 332)]; 'none' -> []."""
    if field == "none":
        return []
    regions = []
    for pair in field.split(","):
        start, end = pair.split("-")
        regions.append((int(start), int(end)))
    return regions


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Collapse overlapping/adjacent (start, end) pairs into contiguous runs.

    post_processing_unitigs.py's non_designable coordinates are raw
    sliding-window hits (one 31bp window per failing position, stepped by
    1bp) -- a long low-complexity stretch can produce thousands of heavily
    overlapping windows for what is really one contiguous bad region.
    Primer3's SEQUENCE_EXCLUDED_REGION has a hard cap on the number of
    intervals it will accept per sequence; passing the raw window list
    through unmerged risks silently overrunning that cap (or just wasting
    thousands of redundant intervals on real ones). Same merge used for the
    masked-bp stats reported alongside POST_PROCESS_MARKERS' results.
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def designable_segments(length: int, non_designable_merged: List[Tuple[int, int]], min_segment_len: int) -> List[Tuple[int, int]]:
    """Complement of the merged non-designable intervals within [0, length),
    dropping any gap shorter than min_segment_len.

    Why segments instead of SEQUENCE_EXCLUDED_REGION: Primer3 caps how many
    excluded intervals it will accept for one template (empirically, real
    7PET markers up to 40kb with ~50% masked coverage blow straight past
    that cap even after merging -- primer3_core rejects them outright with
    "Too many elements for tag SEQUENCE_EXCLUDED_REGION"). Handing primer3 a
    handful of already-clean sub-templates instead sidesteps the cap
    entirely, and is arguably more correct anyway: a 40kb template with half
    of it excluded isn't a meaningful "design somewhere in here" request,
    it's a list of the few places that actually qualify.
    """
    segments = []
    cursor = 0
    for start, end in non_designable_merged:
        if start > cursor:
            segments.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < length:
        segments.append((cursor, length))
    return [(s, e) for s, e in segments if e - s >= min_segment_len]


def header_field(header: str, key: str) -> str:
    """Pull 'key=value' out of a post_processing_unitigs.py FASTA header
    (space-separated 'key=value' tokens after the id)."""
    for token in header.split()[1:]:
        if token.startswith(f"{key}="):
            return token[len(key) + 1 :]
    raise ValueError(f"'{key}=' not found in header: {header!r}")


SEGMENT_ID_SEP = "::"  # marker_id::start-end -- reversed by split_segment_id()


def split_segment_id(sequence_id: str) -> Tuple[str, int]:
    """'2465::15200-18400' -> ('2465', 15200) -- marker id + segment offset,
    so primer3's segment-relative positions can be reported back in the
    original marker's own coordinates."""
    marker_id, _, coords = sequence_id.rpartition(SEGMENT_ID_SEP)
    offset = int(coords.split("-")[0])
    return marker_id, offset


def build_boulder_input(
    fasta_path: Path,
    min_segment_len: int,
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
    """Returns (boulder_text, n_markers, n_segments)."""
    records = []
    n_markers = 0
    n_segments = 0
    for rec in parse_fasta(fasta_path):
        n_markers += 1
        non_designable = merge_intervals(parse_non_designable(header_field(rec.description, "non_designable")))
        for seg_start, seg_end in designable_segments(len(rec.seq), non_designable, min_segment_len):
            n_segments += 1
            segment_seq = rec.seq[seg_start:seg_end].upper()
            block = [
                f"SEQUENCE_ID={rec.id}{SEGMENT_ID_SEP}{seg_start}-{seg_end}",
                f"SEQUENCE_TEMPLATE={segment_seq}",
                "PRIMER_TASK=generic",
                "PRIMER_PICK_LEFT_PRIMER=1",
                "PRIMER_PICK_RIGHT_PRIMER=1",
                "PRIMER_PICK_INTERNAL_OLIGO=0",
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

    return "\n".join(records) + "\n", n_markers, n_segments


##############################################################################
### Boulder-IO output -> TSV


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
    """One row per PRIMER_PAIR_<i>_* pair primer3 actually returned.

    Positions come back relative to the segment primer3 actually saw --
    translated back to the original marker's own coordinates here (+ segment
    offset) so a row is directly usable against the marker FASTA, not a
    sub-template the caller never sees.
    """
    sequence_id = record.get("SEQUENCE_ID", "")
    marker_id, seg_offset = split_segment_id(sequence_id)
    n_returned = int(record.get("PRIMER_PAIR_NUM_RETURNED", "0") or "0")
    rows = []
    for i in range(n_returned):
        left_pos, left_len = record[f"PRIMER_LEFT_{i}"].split(",")
        right_pos, right_len = record[f"PRIMER_RIGHT_{i}"].split(",")
        rows.append(
            {
                "marker_id": marker_id,
                "segment": sequence_id.split(SEGMENT_ID_SEP)[-1],
                "pair_rank": str(i),
                "left_seq": record[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "left_tm": record[f"PRIMER_LEFT_{i}_TM"],
                "left_pos": str(int(left_pos) + seg_offset),
                "left_len": left_len,
                "right_seq": record[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "right_tm": record[f"PRIMER_RIGHT_{i}_TM"],
                "right_pos": str(int(right_pos) + seg_offset),
                "right_len": right_len,
                "product_size": record[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "pair_penalty": record.get(f"PRIMER_PAIR_{i}_PENALTY", ""),
            }
        )
    return rows


TSV_COLUMNS = [
    "marker_id",
    "segment",
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
        fh.write("marker_id\tsegment\treason\n")
        for record in records:
            sequence_id = record.get("SEQUENCE_ID", "")
            marker_id, _ = split_segment_id(sequence_id) if SEGMENT_ID_SEP in sequence_id else (sequence_id, 0)
            segment = sequence_id.split(SEGMENT_ID_SEP)[-1] if SEGMENT_ID_SEP in sequence_id else ""
            reason = record.get("PRIMER_PAIR_EXPLAIN", record.get("PRIMER_ERROR", "no pairs returned"))
            fh.write(f"{marker_id}\t{segment}\t{reason}\n")


##############################################################################
### CLI


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("fasta", type=Path, help="Filtered markers FASTA (POST_PROCESS_MARKERS output)")
    p.add_argument("--label", required=True, help="Prefix for output filenames, e.g. lineage/species id")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--primer3-bin", default="primer3_core", help="primer3_core executable (default: on PATH)")
    p.add_argument(
        "--min-segment-length",
        type=int,
        default=100,
        help="Shortest contiguous non-masked stretch worth submitting to primer3 (default: 100bp)",
    )
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

    boulder_input, n_markers, n_segments = build_boulder_input(
        args.fasta,
        min_segment_len=args.min_segment_length,
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
    if n_segments == 0:
        sys.exit(
            f"{n_markers:,} marker(s) read, but none had a contiguous non-masked stretch "
            f">= {args.min_segment_length}bp -- nothing to submit to primer3."
        )

    print(
        f"Running primer3_core on {n_segments:,} clean segment(s) from {n_markers:,} marker(s) "
        f"in {args.fasta} ...",
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
            f"(a single unparseable segment breaks the whole Boulder-IO batch):\n{result.stderr}"
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
    print(f"  Markers submitted          : {n_markers:,}", file=sys.stderr)
    print(f"  Clean segments submitted   : {n_segments:,}", file=sys.stderr)
    print(f"  Markers with >=1 pair      : {n_with_primers:,} ({100 * n_with_primers / n_markers:.1f}%)", file=sys.stderr)
    print(f"  Markers with 0 pairs       : {n_markers - n_with_primers:,}", file=sys.stderr)
    print(f"  Total primer pairs found   : {len(all_rows):,}", file=sys.stderr)
    print(f"  Wrote {primers_path}", file=sys.stderr)
    print(f"  Wrote {no_primers_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
