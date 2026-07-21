#!/usr/bin/env python3

"""
LINEAGE-CORE / CATCH-ALL FILTERING PIPELINE
(species-agnostic, lineage-agnostic, background-db-agnostic)

Builds a lineage-specific marker k-mer set (unitigs) for a target lineage
within a target species, refined so the markers are specific both within
the species (not shared with other lineages) and against the wider genomic
background (not shared with large external reference collections such as
ATB or GTDB).


TERMINOLOGY
  species-agnostic       A step's logic doesn't reference any specific
                          species by name — the target species is entirely
                          defined by the SPECIES_INDEX parameter passed in.
  lineage-agnostic        Same, for lineage — no lineage ID is hardcoded,
                          all lineage specifics come from LINEAGE_LIST /
                          LINEAGE_MAP.
  background-db-agnostic  The pipeline doesn't care which reference database
                          (ATB, GTDB, or any future db) is used as
                          background — any db can be plugged in via
                          BACKGROUND_DBS without code changes.
  index-native            The operation is performed directly by the SBWT /
                          Themisto2 CLI tools on a real k-mer index (e.g.
                          via `sbwt difference`), not derived in Python.
  Python-derived          The object was built in this script via boolean-
                          vector arithmetic over colour-set membership
                          (cheap, avoids building an intermediate index),
                          and is therefore NOT a real SBWT/Themisto2 index
                          until it is explicitly rebuilt as one.


Parameters (passed in, never hardcoded):
  SPECIES_INDEX       = combined SBWT/Themisto index for target species
                        (e.g. S. pneumoniae combined index, unitigs + colors)
  LINEAGE_LIST        = list of lineage IDs to process (e.g. all GPSCs of interest)
  LINEAGE_MAP         = mapping of lineage ID -> list of genome/color IDs
  BACKGROUND_DBS      = list of background reference indexes to check against
                        (e.g. [ATB, GTDB] — extensible to any future db)
  MIN_FREQ            = minimum presence fraction a unitig must reach within
                        the lineage to be kept. This is the one real
                        parameter; CORE_MODE / RELAXED_MODE / CATCHALL_MODE
                        are just named presets (profiles) for it — e.g.
                        ~0.95 for "core", ~0.0-0.1 for "catch-all" — always
                        directly overridable regardless of preset chosen.
  SPECIFICITY_CUTOFF  = max allowed presence outside lineage (e.g. 0.1)


DEFINITIONS (A -> G, in the order they appear in the pipeline)

  A        = species-wide SBWT/Themisto2 index — all genomes of the target
             species, colored. Input, built once, reused for every lineage.

  B        = per-lineage SBWT index — the subset of A restricted to one
             lineage's genomes. Input; only built if doing an index-native
             D (Step 2) — otherwise B is just the lineage's colour IDs
             within A, no separate index required.

  ATB/GTDB = background reference SBWT indexes, built independently from
             large external genome collections (Themisto2-exported
             unitigs). Not derived from A or B.

  C        = ATB − A   (and C_gtdb = GTDB − A)
             "background exclusion set" — k-mers/unitigs that exist in the
             wider world but are NOT part of the target species at all.
             Species-agnostic: depends only on A, not on any lineage, so
             it's computed once per run and reused across every lineage.

  D        = A − B
             "cross-lineage background" — k-mers/unitigs found somewhere
             else in the species but NOT in this lineage. Per-lineage,
             index-native.

  E        = min_freq-filtered candidates from B
             "lineage-core / catch-all candidates" — unitigs kept because
             their presence across the lineage's own genomes crosses
             MIN_FREQ (mode supplies the default cutoff — ~0.95 for "core",
             ~0.0-0.1 for "catch-all" — always overridable via --min-freq).
             Built in Python from a boolean colour-membership vector, per
             lineage — Python-derived, not yet a real index. THIS SCRIPT'S
             OUTPUT STOPS HERE.

  F        = E − D   (computed after rebuilding E as a real index — Step 4a)
             "candidate lineage-specific k-mers" — E's candidates with
             anything that also occurs in other lineages of the same
             species (D) removed. Confirms specificity WITHIN the species.

  G        = F − C   (and G_gtdb = F − C_gtdb)
             "candidate lineage-specific AND background-specific k-mers" —
             F's candidates with anything that also occurs in the wider
             background database (ATB / GTDB) removed. Confirms
             specificity against the outside world, not just within the
             species. THIS IS THE FINAL OUTPUT: the lineage marker set.


PIPELINE — four steps, in order

  Step 1 — Background exclusion (species-agnostic, computed once)
    C      = ATB  − A     [sbwt difference ATB.sbwt  A.sbwt -o C.sbwt]
    C_gtdb = GTDB − A     [sbwt difference GTDB.sbwt A.sbwt -o C_gtdb.sbwt]
    Runs once per pipeline run, independent of lineage — reused for every
    lineage processed afterwards.

  Step 2 — Cross-lineage background (per lineage, index-native)
    D = A − B             [sbwt difference A.sbwt B.sbwt -o D.sbwt]
    Only needed if doing an index-native D for benchmarking/validation
    (see GPSC1 benchmark ticket). A Python-derived version of D can be
    computed from the same boolean-vector logic used to build E (Step 3),
    but that is benchmark-only, valid only at the threshold=100% boundary.

  Step 3 — Lineage-core / catch-all min_freq candidates (per lineage,
           Python, boolean vector — THIS SCRIPT'S SCOPE)
    E = min_freq-filtered B unitigs
        (mode presets the cutoff — core: ~0.95; catch-all: ~0.0-0.1 — via
         sum(a & lineage_membership) / sum(lineage_membership); always
         overridable directly via --min-freq)
    This script (core_catchall_filter.py) writes E's surviving unitigs to
    a FASTA and stops. Everything in Step 4 happens outside this file.

  Step 4 — Candidate refinement (rebuild + index-native diff, external)
    4a. Rebuild E's unitig FASTA into a real SBWT + Themisto2 index
        (E.sbwt). Necessary because E is Python-derived (Step 3), not a
        native index — it has to become one before it can be diffed
        index-natively.
    4b. F      = E − D      [sbwt difference E.sbwt D.sbwt -o F.sbwt]
    4c. G      = F − C      [sbwt difference F.sbwt C.sbwt -o G.sbwt]
        G_gtdb = F − C_gtdb [sbwt difference F.sbwt C_gtdb.sbwt -o G_gtdb.sbwt]


OPEN QUESTION
  Whether G_gtdb replaces, supplements, or runs parallel to G (the
  ATB-only result) is not yet decided — i.e. whether GTDB should be
  checked instead of, in addition to, or independently from ATB. This is
  a scope question, not an implementation detail.


NOTES ON GENERALISATION
  - Nothing above references "S. pneumoniae," "GPSC," "ATB," or "GTDB" by
    name — all species/lineage/db specifics live in SPECIES_INDEX,
    LINEAGE_LIST, LINEAGE_MAP, and BACKGROUND_DBS, passed as
    parameters/config.
  - Adding a new background db (e.g. a third reference beyond ATB/GTDB)
    only requires adding it to BACKGROUND_DBS — Steps 1 and 4c loop over
    it automatically, no code changes needed.
  - Adding a new target species requires only a new SPECIES_INDEX and its
    own LINEAGE_LIST/LINEAGE_MAP — same script, same logic.
  - `sbwt difference` does the heavy, index-native subtraction (C, D, F, G)
    — this script does not reimplement that. It only builds E.


NOTE ON FILE SCOPE / CURRENT STATUS
  This file (core_catchall_filter.py) contains only the Step 3 logic:
  computing per-lineage core/outside presence and writing the
  threshold-filtered candidate set E. Steps 1, 2, and 4 (all native
  `sbwt difference` calls, plus the Step 4a index rebuild) are not
  implemented here and are not planned to be — they run as separate
  CLI/sbwt/Themisto2 commands outside this script.
  The existing species-exclusion script (used for PAT-3461, producing C
  via colour-based filtering) stays separate, unchanged, in
  filter_species_from_atb.py.
"""

import argparse
import gzip
import sys
from pathlib import Path
import numpy as np

##############################################################################
### I/O helper

def _open(path):
    path = Path(path)
    if path.exists():
        return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)
    gz_path = Path(str(path) + ".gz")
    if gz_path.exists():
        return gzip.open(gz_path, "rt")
    no_gz_path = Path(str(path)[:-3]) if str(path).endswith(".gz") else None
    if no_gz_path and no_gz_path.exists():
        return open(no_gz_path)
    raise FileNotFoundError(f"Could not find file: {path} (also tried .gz variant)")


def _open_out(path):
    """Open a file for writing, gzip-compressed if the path ends in .gz."""
    path = Path(path)
    return gzip.open(path, "wt") if path.suffix == ".gz" else open(path, "w")

##############################################################################
### I/O loading

def load_unitig_colorset_ids(index_export_path) -> np.ndarray:
    colorset_ids = []

    with _open(index_export_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                fields = dict(
                    field.split("=") for field in header.split()
                )
                colorset_ids.append(int(fields["color_set_id"]))

    return np.array(colorset_ids, dtype=np.int64)

def load_lineage_map(lineage_mapping_path) -> dict[str, list[int]]:
    """
    Parse a TSV with columns Sample_ID, lineage_id into a
    lineage_id -> [color_ids] dict.

    color_id is 0-indexed and equals the row's position in the file
    (header on line 1, data starts line 2), matching the build-time
    genome order used for color ID assignment.
    """
    lineage_map: dict[str, list[int]] = {}
    with _open(lineage_mapping_path) as fh:
        next(fh)  # skip header row
        for color_id, line in enumerate(fh):
            sample_id, lineage_id = line.rstrip("\n").split("\t")
            lineage_map.setdefault(lineage_id, []).append(color_id)
    return lineage_map

##############################################################################
### Per-lineage membership and fraction calculations:

def build_lineage_color_mask(lineage_color_ids, n_total_colors) -> np.ndarray:
    """
    Build a boolean vector of length n_total_colors, True at positions
    corresponding to genomes belonging to the target lineage.
    """
    lineage_mask = np.zeros(n_total_colors, dtype=np.bool_)
    lineage_mask[lineage_color_ids] = True
    return lineage_mask

def compute_colorset_fractions_streaming(
    color_sets_path,
    lineage_masks: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """
    Single streaming pass over export.color_sets.txt(.gz) that computes, for
    every requested lineage at once, both the within-lineage and
    outside-lineage presence fraction of each color_set_id:
        within  = |color_set ∩ lineage|  / |lineage|
        outside = |color_set ∩ ~lineage| / |~lineage|

    outside_pct is what SPECIFICITY_CUTOFF (see module docstring) refers
    to — how common this k-mer is *outside* the target lineage. Computed
    here (not just within_pct) so callers can score candidates by how
    cleanly they separate the lineage from everything else, not just
    threshold on within-lineage presence alone.

    Unlike parsing the file into a color_set_id -> [color_ids] dict up
    front, this never retains the per-color_set color_ids beyond the line
    they came from — only the resulting per-lineage fractions are kept.
    Peak memory is O(n_colorsets * n_lineages) instead of O(total color_ids
    in the file).

    Returns
    -------
    tuple[dict[str, np.ndarray], dict[str, np.ndarray]]
        (within_fractions_by_lineage, outside_fractions_by_lineage), each
        lineage_id -> dense array of presence fractions indexed by
        color_set_id.
    """
    lineage_sizes = {
        lineage_id: int(np.sum(mask)) for lineage_id, mask in lineage_masks.items()
    }
    for lineage_id, size in lineage_sizes.items():
        if size == 0:
            raise ValueError(f"lineage '{lineage_id}' mask has zero members — empty lineage?")

    n_total_colors = next(iter(lineage_masks.values())).shape[0]
    outside_masks = {lineage_id: ~mask for lineage_id, mask in lineage_masks.items()}
    outside_sizes = {
        lineage_id: n_total_colors - size for lineage_id, size in lineage_sizes.items()
    }
    for lineage_id, size in outside_sizes.items():
        if size == 0:
            raise ValueError(f"lineage '{lineage_id}' has no genomes outside it — "
                              f"cannot compute an outside-lineage fraction")

    seen_colorset_ids = []
    within_by_lineage = {lineage_id: [] for lineage_id in lineage_masks}
    outside_by_lineage = {lineage_id: [] for lineage_id in lineage_masks}

    with _open(color_sets_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            color_set_id = int(parts[0].split("=")[1])
            color_ids = np.array(
                [int(tok) for tok in parts[2:]], dtype=np.int64
            )  # skip color_set_id=N, size=M; discarded after this iteration

            seen_colorset_ids.append(color_set_id)
            for lineage_id, mask in lineage_masks.items():
                within_by_lineage[lineage_id].append(
                    np.sum(mask[color_ids]) / lineage_sizes[lineage_id]
                )
                outside_by_lineage[lineage_id].append(
                    np.sum(outside_masks[lineage_id][color_ids]) / outside_sizes[lineage_id]
                )

    colorset_ids_arr = np.array(seen_colorset_ids, dtype=np.int64)
    max_colorset_id = colorset_ids_arr.max() if len(colorset_ids_arr) else -1

    def _to_dense(by_lineage):
        dense_by_lineage = {}
        for lineage_id, values in by_lineage.items():
            dense = np.zeros(max_colorset_id + 1, dtype=np.float64)
            dense[colorset_ids_arr] = np.array(values, dtype=np.float64)
            dense_by_lineage[lineage_id] = dense
        return dense_by_lineage

    return _to_dense(within_by_lineage), _to_dense(outside_by_lineage)

def map_unitig_fractions(
    unitig_colorset_ids: np.ndarray,
    fraction_array: np.ndarray,
) -> np.ndarray:
    return fraction_array[unitig_colorset_ids]

##############################################################################
### Mode presets
#
# mode is a named preset (profile) for the one real parameter, min_freq —
# the minimum presence fraction a unitig must reach within the lineage to
# be kept. Each mode just supplies min_freq's default; --min-freq always
# overrides it directly, regardless of mode.

COREMODE_MIN_FREQ = 0.95
RELAXMODE_MIN_FREQ = 0.5
CATCHMODE_MIN_FREQ = 0.0

##############################################################################
### Thresholding logic

def min_freq_mask(
    fractions: np.ndarray,
    mode: str,
    min_freq: float | None = None,
) -> np.ndarray:
    """
    mode='core':     fractions >= min_freq   (default min_freq = COREMODE_MIN_FREQ)
                      STRICT — k-mer must be present in ~95%+ of the
                      lineage's genomes. Near-universal, tolerant of a
                      small number of dropouts (assembly gaps,
                      sequencing noise, one-off bad samples).

    mode='relaxed':   fractions >  min_freq   (default min_freq = RELAXMODE_MIN_FREQ)
                      MAJORITY — k-mer must be present in >50% of the
                      lineage's genomes. Includes core + common
                      accessory content, excludes rare/sporadic k-mers.

    mode='catchall':  fractions >  min_freq   (default min_freq = CATCHMODE_MIN_FREQ)
                      ANY — k-mer must be present in at least one
                      genome of the lineage. Core + all accessory. The
                      strict '>' means a unitig absent from every genome
                      in the lineage (fraction == 0) never passes, even
                      at the default floor of 0.0.

    min_freq overrides the mode's default cutoff.
    """
    if mode == "core":
        cutoff = min_freq if min_freq is not None else COREMODE_MIN_FREQ
        return fractions >= cutoff
    elif mode == "relaxed":
        cutoff = min_freq if min_freq is not None else RELAXMODE_MIN_FREQ
        return fractions > cutoff
    elif mode == "catchall":
        cutoff = min_freq if min_freq is not None else CATCHMODE_MIN_FREQ
        return fractions > cutoff
    else:
        raise ValueError(f"unknown mode: {mode!r}")

##############################################################################
###  Output writing

def write_lineage_unitig_ids(mask: np.ndarray, out_path):
    """Write unitig_ids where mask is True to out_path, one per line."""
    passing_ids = np.nonzero(mask)[0]
    total = len(mask)
    kept_count = len(passing_ids)
    dropped_count = total - kept_count

    with _open_out(out_path) as out_fh:
        for unitig_id in passing_ids:
            out_fh.write(f"{unitig_id}\n")

    return total, kept_count, dropped_count

def write_lineage_summary_stats(lineage_id, mode, min_freq, total, kept_count,
                                 dropped_count, stats_out_path):
    """Write a summary of one lineage+mode run to stats_out_path."""
    stats_lines = [
        f"{'Lineage':<21} : {lineage_id}",
        f"{'Mode':<21} : {mode}",
        f"{'Min freq':<21} : {min_freq}",
        f"{'Total unitigs scanned':<21} : {total:,}",
        f"{'Kept (E)':<21} : {kept_count:,}",
        f"{'Dropped':<21} : {dropped_count:,}",
    ]
    Path(stats_out_path).write_text("\n".join(stats_lines) + "\n")

def write_specificity_scores(mask, core_fractions, outside_fractions, out_path):
    """
    Write per-unitig specificity scores for the kept (E) set, i.e. unitigs
    where mask is True.

    score = core_pct - outside_pct

    Mathematically identical to Youden's J statistic
    (sensitivity + specificity - 1), the standard diagnostic-marker-quality
    metric: +1 = perfect marker (present in 100% of the target lineage, 0%
    of everything else); 0 = no discriminating power; negative = the k-mer
    is actually more common outside the lineage than inside it.
    """
    passing_ids = np.nonzero(mask)[0]
    with _open_out(out_path) as out_fh:
        out_fh.write("unitig_id\tcore_pct\toutside_pct\tspecificity_score\n")
        for unitig_id in passing_ids:
            core_pct = core_fractions[unitig_id]
            outside_pct = outside_fractions[unitig_id]
            out_fh.write(
                f"{unitig_id}\t{core_pct:.6f}\t{outside_pct:.6f}\t{core_pct - outside_pct:.6f}\n"
            )

##############################################################################
### Orchestrate:

def process_lineage(
    lineage_id: str,
    fraction_array: np.ndarray,
    unitig_colorset_ids: np.ndarray,
    mode: str,
    min_freq: float | None,
    out_path,
    stats_out_path=None,
    outside_fraction_array: np.ndarray | None = None,
    specificity_out_path=None,
):
    unitig_fractions = map_unitig_fractions(unitig_colorset_ids, fraction_array)

    mask = min_freq_mask(unitig_fractions, mode, min_freq)

    total, kept_count, dropped_count = write_lineage_unitig_ids(mask, out_path)

    if stats_out_path:
        write_lineage_summary_stats(
            lineage_id, mode, min_freq, total, kept_count, dropped_count,
            stats_out_path,
        )

    if specificity_out_path:
        if outside_fraction_array is None:
            raise ValueError("specificity_out_path given but outside_fraction_array is None")
        unitig_outside_fractions = map_unitig_fractions(unitig_colorset_ids, outside_fraction_array)
        write_specificity_scores(mask, unitig_fractions, unitig_outside_fractions, specificity_out_path)

    return total, kept_count, dropped_count

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-lineage core presence for unitigs in a "
                     "Themisto2-exported species index, then apply a "
                     "threshold to build a filtered candidate marker set "
                     "(E), written as unitig_id lists.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--unitigs", required=True, type=Path,
                         help="export.unitigs.fa(.gz) from Themisto2 export "
                              "(species-wide, colored). Used to derive "
                              "unitig_colorset_ids.")
    parser.add_argument("--color-sets", required=True, type=Path,
                         help="export.color_sets.txt(.gz) from Themisto2 export. "
                              "Format per line: "
                              "'color_set_id=N size=M c1 c2 c3 ...'")
    parser.add_argument("--lineage-mapping", required=True, type=Path,
                         help="TSV(.gz) with columns Sample_ID, lineage_id, in "
                              "the same genome order used at index build time "
                              "(color_id = row index, 0-indexed, header on "
                              "line 1).")
    parser.add_argument("--n-total-colors", required=True, type=int,
                         help="Total number of colors/genomes in the species "
                              "index (from export.metadata.txt's num_colors).")
    parser.add_argument("--lineages", required=True, nargs="+",
                         help="One or more target lineage IDs to filter for, "
                              "e.g. --lineages 1 2 3 (GPSC1, GPSC2, GPSC3). "
                              "Must match lineage_id values in "
                              "--lineage-mapping.")
    parser.add_argument("--out-dir", required=True, type=Path,
                         help="Directory to write per-lineage unitig_id list "
                              "files (E). One file per lineage, named "
                              "'{lineage_id}_{mode}_unitig_ids.txt', plain "
                              "text with one unitig_id per line. Created if "
                              "it doesn't exist.")
    parser.add_argument("--stats-out-dir", type=Path, default=None,
                         help="Optional directory to write per-lineage "
                              "summary stats (total/kept/dropped counts), "
                              "named the same way as --out-dir. Created if "
                              "it doesn't exist. If omitted, no stats are "
                              "written.")
    parser.add_argument("--specificity-out-dir", type=Path, default=None,
                         help="Optional directory to write per-unitig "
                              "specificity scores for the kept (E) set: "
                              "core_pct, outside_pct, and "
                              "specificity_score = core_pct - outside_pct "
                              "(Youden's J -- diagnostic marker quality; "
                              "+1 = perfect, 0 = no discriminating power, "
                              "negative = more common outside the lineage "
                              "than inside it). Named "
                              "'{lineage_id}_{mode}_specificity.tsv'. "
                              "Created if it doesn't exist. If omitted, not "
                              "written.")
    parser.add_argument("-m", "--mode",
                         choices=["core", "relaxed", "catchall"],
                         default="core",
                         help="Preset for --min-freq. 'core': strict, "
                              f"present in ~95%% of lineage genomes "
                              f"(default min_freq {COREMODE_MIN_FREQ}). "
                              "'relaxed': majority, present in >50%% of "
                              f"lineage genomes (default min_freq "
                              f"{RELAXMODE_MIN_FREQ}). 'catchall': any, "
                              "present in at least one lineage genome "
                              f"(default min_freq {CATCHMODE_MIN_FREQ}). "
                              "All overridable with --min-freq.")
    parser.add_argument("-f", "--min-freq",
                         type=float,
                         default=None,
                         help="Minimum presence fraction a unitig must "
                              "reach within the lineage to be kept. "
                              "Overrides the mode's default "
                              f"({COREMODE_MIN_FREQ} core / "
                              f"{RELAXMODE_MIN_FREQ} relaxed / "
                              f"{CATCHMODE_MIN_FREQ} catchall). Must be "
                              "between 0.0 and 1.0.")
    return parser.parse_args()

def validate_min_freq_args(args):
    if args.min_freq is not None and not (0.0 <= args.min_freq <= 1.0):
        raise ValueError(f"--min-freq must be between 0.0 and 1.0, got {args.min_freq}")

def main():
    args = parse_args()
    validate_min_freq_args(args)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.stats_out_dir:
        args.stats_out_dir.mkdir(parents=True, exist_ok=True)
    if args.specificity_out_dir:
        args.specificity_out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading unitig colorset IDs from {args.unitigs} ...", file=sys.stderr)
    unitig_colorset_ids = load_unitig_colorset_ids(args.unitigs)
    print(f"  {len(unitig_colorset_ids):,} unitig(s) loaded", file=sys.stderr)

    print(f"Loading lineage mapping from {args.lineage_mapping} ...", file=sys.stderr)
    lineage_map = load_lineage_map(args.lineage_mapping)
    print(f"  {len(lineage_map):,} lineage(s) found in mapping", file=sys.stderr)

    lineage_masks = {}
    for lineage_id in args.lineages:
        if lineage_id not in lineage_map:
            print(f"WARNING: lineage '{lineage_id}' not found in mapping, skipping",
                  file=sys.stderr)
            continue
        lineage_masks[lineage_id] = build_lineage_color_mask(
            lineage_map[lineage_id], args.n_total_colors
        )

    print(f"Streaming color sets from {args.color_sets} ...", file=sys.stderr)
    fractions_by_lineage, outside_fractions_by_lineage = compute_colorset_fractions_streaming(
        args.color_sets, lineage_masks
    )
    print(f"  presence fractions computed for {len(fractions_by_lineage):,} lineage(s)",
          file=sys.stderr)

    for lineage_id in lineage_masks:
        out_path = args.out_dir / f"{lineage_id}_{args.mode}_unitig_ids.txt"
        stats_out_path = (
            args.stats_out_dir / f"{lineage_id}_{args.mode}_unitig_ids.txt"
            if args.stats_out_dir else None
        )
        specificity_out_path = (
            args.specificity_out_dir / f"{lineage_id}_{args.mode}_specificity.tsv"
            if args.specificity_out_dir else None
        )

        print(f"Processing lineage '{lineage_id}' (mode={args.mode}) ...",
              file=sys.stderr)
        total, kept_count, dropped_count = process_lineage(
            lineage_id=lineage_id,
            fraction_array=fractions_by_lineage[lineage_id],
            unitig_colorset_ids=unitig_colorset_ids,
            mode=args.mode,
            min_freq=args.min_freq,
            out_path=out_path,
            stats_out_path=stats_out_path,
            outside_fraction_array=outside_fractions_by_lineage[lineage_id],
            specificity_out_path=specificity_out_path,
        )
        print(f"  Total: {total:,}  Kept (E): {kept_count:,}  Dropped: {dropped_count:,}",
              file=sys.stderr)


if __name__ == "__main__":
    main()