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
  CORE_THRESHOLD      = core presence cutoff (e.g. 0.9 for "core", 0.0-0.1 for "catch-all")
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

  E        = threshold-filtered candidates from B
             "lineage-core / catch-all candidates" — unitigs kept because
             their presence across the lineage's own genomes crosses
             CORE_THRESHOLD (~0.9-1.0 for "core" mode, ~0.0-0.1 for
             "catch-all" mode). Built in Python from a boolean colour-
             membership vector, per lineage — Python-derived, not yet a
             real index. THIS SCRIPT'S OUTPUT STOPS HERE.

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

  Step 3 — Lineage-core / catch-all threshold candidates (per lineage,
           Python, boolean vector — THIS SCRIPT'S SCOPE)
    E = threshold-filtered B unitigs
        (core mode: CORE_THRESHOLD ~0.9-1.0; catch-all mode: ~0.0-0.1;
         via sum(a & lineage_membership) / sum(lineage_membership))
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
import numpy

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

def build_colorset_to_colorids():
    pass

def load_lineage_map():
    pass
##############################################################################
### Per-lineage membership and fraction calculations:
def build_lineage_color_mask():
    pass
def compute_colorset_presence_fractions():
    pass
def map_unitig_fractions():
    pass

##############################################################################
### Thresholding logic
def threshold_mask(
    fractions: np.ndarray,
    mode: str,
    threshold: float | None = None,
) -> np.ndarray:
    """
    mode='core':     fractions >= threshold   (default threshold = 0.95)
                      STRICT — k-mer must be present in ~95%+ of the
                      lineage's genomes. Near-universal, tolerant of a
                      small number of dropouts (assembly gaps,
                      sequencing noise, one-off bad samples).

    mode='relaxed':   fractions >  threshold   (default threshold = 0.5)
                      MAJORITY — k-mer must be present in >50% of the
                      lineage's genomes. Includes core + common
                      accessory content, excludes rare/sporadic k-mers.

    mode='catchall':  fractions >  threshold   (default threshold = 0.0)
                      ANY — k-mer must be present in at least one
                      genome of the lineage. Core + all accessory,
                      no noise floor.

    threshold overrides the default cutoff for any mode.
    """
    if mode == "core":
        cutoff = threshold if threshold is not None else 0.95
        return fractions >= cutoff
    elif mode == "relaxed":
        cutoff = threshold if threshold is not None else 0.5
        return fractions > cutoff
    elif mode == "catchall":
        cutoff = threshold if threshold is not None else 0.0
        return fractions > cutoff
    else:
        raise ValueError(f"unknown mode: {mode!r}")

def validate_threshold_args(args):
    if args.threshold is not None and not (0.0 <= args.threshold <= 1.0):
        raise ValueError(f"--threshold must be between 0.0 and 1.0, got {args.threshold}")
##############################################################################
###  Output writing
def write_lineage_summary_stats():
    pass

##############################################################################
### Orchestrate:
def process_lineage():
    pass

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-lineage core presence and outside-lineage "
                     "specificity for unitigs in a Themisto2-exported species "
                     "index, then apply a threshold to build a filtered "
                     "candidate marker set (E).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--unitigs", required=True, type=Path,
                         help="export.unitigs.fa from Themisto2 export "
                              "(species-wide, colored).")
    parser.add_argument("--color-sets", required=True, type=Path,
                         help="export.color_sets.txt(.gz) from Themisto2 export.")
    parser.add_argument("--lineage-mapping", required=True, type=Path,
                         help="TSV with columns Sample_ID, lineage_id, in the "
                              "same genome order used at index build time "
                              "(colour_id = row index, 0-indexed).")
    parser.add_argument("--n-total-colors", required=True, type=int,
                         help="Total number of colours/genomes in the species "
                              "index (from export.metadata.txt's num_colors).")
    parser.add_argument("--lineage", required=True,
                         help="Target lineage ID to filter for, e.g. '1' (GPSC1).")
    parser.add_argument("--core-threshold", required=True, type=float,
                         help="Minimum core_pct required to keep a unitig. "
                              "Use ~0.9-1.0 for core mode, ~0.0-0.1 for "
                              "catch-all mode.")
    parser.add_argument("--specificity-cutoff", type=float, default=None,
                         help="Optional maximum outside_pct allowed. If not "
                              "given, specificity filtering is skipped at "
                              "this stage.")
    parser.add_argument("--out", required=True, type=Path,
                         help="Output path for filtered unitig FASTA (E). "
                              "Use .gz suffix to compress.")
    parser.add_argument("--stats-out", type=Path, default=None,
                         help="Optional path to write summary stats "
                              "(total/kept/dropped counts).")
    parser.add_argument("-m", "--mode", 
                         choices=["core", "catchall"],
                         default="core",
                         help="'core': strict, present in 95%% of lineage genomes (override with --threshold). "
                              "'catchall': relaxed, present in >50%% of lineage genomes (override with --threshold).",)
    parser.add_argument("-t", "--threshold", 
                         type=float, 
                         default=None,
                         help="Override default cutoff (0.95 core / 0.5 catchall).")
    return parser.parse_args()

def validate_threshold_args():
    pass
def main():
    pass

##############################################################################
##############################################################################
##############################################################################
### Old version to base new functions on
def build_colourid_to_lineage(lineage_mapping_path):
    """
    Build a mapping of colour_id -> lineage_id from a TSV file with columns:
      Sample_ID, lineage_id

    Colour IDs are always 0-indexed (Themisto2/SBWT convention), matching
    the row order of `lineage_mapping_path`, which has been confirmed to
    mirror the build-time genome order used for colour ID assignment.
    (header on line 1, data starts line 2 -> colour_id = row_index)
    """
    colourid_to_lineage = {}
    with open(lineage_mapping_path) as fh:
        next(fh)  # skip header row
        for colour_id, line in enumerate(fh):
            sample_id, lineage = line.rstrip("\n").split("\t")
            colourid_to_lineage[colour_id] = lineage
    return colourid_to_lineage


def find_lineage_colour_ids(colourid_to_lineage, lineage):
    """
    Return the set of colour IDs assigned to `lineage`.
    """
    matched_ids = {cid for cid, lin in colourid_to_lineage.items() if lin == lineage}
    if not matched_ids:
        raise ValueError(f"No colour IDs matched lineage '{lineage}'")
    return matched_ids

def build_colorset_to_colourids(color_sets_path):
    """
    Parse export.color_sets.txt(.gz) into a color_set_id -> [colour_ids] dict.

    Format per line: "color_set_id=N size=M c1 c2 c3 ..."
    """
    colorset_to_colourids = {}
    with _open(color_sets_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            color_set_id = int(parts[0].split("=")[1])
            colour_ids = [int(tok) for tok in parts[2:]]  # skip color_set_id=N, size=M
            colorset_to_colourids[color_set_id] = colour_ids
    return colorset_to_colourids

def build_lineage_membership_vector(lineage_colour_ids, n_total_colors):
    """
    Build a boolean vector of length n_total_colors, True at positions
    corresponding to genomes belonging to the target lineage.
    """
    lineage_membership = numpy.zeros(n_total_colors, dtype=numpy.bool_)
    lineage_membership[list(lineage_colour_ids)] = True
    return lineage_membership

def iter_unitig_lineage_stats(unitigs_path, colorset_to_colourids, lineage_membership):
    """
    Stream export.unitigs.fa one record at a time. For each unitig, look up
    its colour set (via its color_set_id) and compute:
      - core_pct    : fraction of lineage genomes where this unitig is present
      - outside_pct : fraction of non-lineage genomes where this unitig is
                       also present (specificity / leakage)

    Yields (header, seq, color_set_id, core_pct, outside_pct) per unitig.

    n_total_colors is inferred from lineage_membership's length, since it
    must match anyway (both index into the same colour space).
    """
    n_total_colors = len(lineage_membership)
    lineage_size = numpy.sum(lineage_membership)
    outside_size = numpy.sum(numpy.invert(lineage_membership))

    with _open(unitigs_path) as in_fh:
        current_header = None
        current_colorset = None
        current_seq = None

        def build_stats():
            colour_ids = colorset_to_colourids[current_colorset]
            a = numpy.zeros(n_total_colors, dtype=numpy.bool_)
            a[colour_ids] = True

            core_pct = numpy.sum(a & lineage_membership) / lineage_size
            outside_pct = numpy.sum(a & numpy.invert(lineage_membership)) / outside_size

            return current_header, current_seq, current_colorset, core_pct, outside_pct

        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    yield build_stats()
                current_header = line[1:].strip()
                current_colorset = None
                for tok in current_header.split():
                    if tok.startswith("color_set_id="):
                        current_colorset = int(tok.split("=", 1)[1])
                        break
                current_seq = None
            else:
                current_seq = line

        if current_header is not None:
            yield build_stats()  # last record



def write_threshold_filtered_fasta(unitigs_path, colorset_to_colourids,
                                    lineage_membership, core_threshold,
                                    specificity_cutoff, out_path,
                                    stats_out_path=None):
    """
    Stream export.unitigs.fa, applying threshold filtering to build E:
    keep a unitig if its core_pct meets core_threshold, and (optionally)
    its outside_pct stays below specificity_cutoff.

    core_threshold semantics:
      - core mode:      core_threshold ~ 0.9-1.0 (present in most/all of lineage)
      - catch-all mode: core_threshold ~ 0.0-0.1 (present in any/some of lineage)
    In both modes the same comparison (core_pct >= core_threshold) applies —
    the caller decides which mode by choosing the threshold value.

    specificity_cutoff is optional (pass None to skip it, keeping this
    function usable for both "core only" and "core + specificity" filtering
    without a separate code path).

    Returns (total, kept_count, dropped_count).
    """
    total = 0
    kept_count = 0
    dropped_count = 0

    with _open_out(out_path) as out_fh:
        for header, seq, color_set_id, core_pct, outside_pct in iter_unitig_lineage_stats(
            unitigs_path, colorset_to_colourids, lineage_membership
        ):
            total += 1

            passes_core = core_pct >= core_threshold
            passes_specificity = (
                specificity_cutoff is None or outside_pct < specificity_cutoff
            )

            if passes_core and passes_specificity:
                out_fh.write(f">{header}\n{seq}\n")
                kept_count += 1
            else:
                dropped_count += 1

    if stats_out_path:
        stats_lines = [
            f"Total unitigs scanned : {total:,}",
            f"Kept (E)              : {kept_count:,}",
            f"Dropped               : {dropped_count:,}",
            f"Core threshold        : {core_threshold}",
            f"Specificity cutoff    : {specificity_cutoff}",
        ]
        Path(stats_out_path).write_text("\n".join(stats_lines) + "\n")

    return total, kept_count, dropped_count


def main():
    args = parse_args()

    print(f"Loading lineage mapping from {args.lineage_mapping} ...", file=sys.stderr)
    colourid_to_lineage = build_colourid_to_lineage(args.lineage_mapping)

    print(f"Finding colour IDs for lineage '{args.lineage}' ...", file=sys.stderr)
    lineage_colour_ids = find_lineage_colour_ids(colourid_to_lineage, args.lineage)
    print(f"  {len(lineage_colour_ids):,} genome(s) in lineage '{args.lineage}'",
          file=sys.stderr)

    print("Building lineage membership vector ...", file=sys.stderr)
    lineage_membership = build_lineage_membership_vector(
        lineage_colour_ids, args.n_total_colors
    )

    print(f"Parsing colour sets from {args.color_sets} ...", file=sys.stderr)
    colorset_to_colourids = build_colorset_to_colourids(args.color_sets)
    print(f"  {len(colorset_to_colourids):,} colour set(s) loaded", file=sys.stderr)

    print(f"Filtering unitigs (core_threshold={args.core_threshold}, "
          f"specificity_cutoff={args.specificity_cutoff}) ...", file=sys.stderr)
    total, kept_count, dropped_count = write_threshold_filtered_fasta(
        unitigs_path=args.unitigs,
        colorset_to_colourids=colorset_to_colourids,
        lineage_membership=lineage_membership,
        core_threshold=args.core_threshold,
        specificity_cutoff=args.specificity_cutoff,
        out_path=args.out,
        stats_out_path=args.stats_out,
    )

    print(f"Done. Total: {total:,}  Kept (E): {kept_count:,}  "
          f"Dropped: {dropped_count:,}", file=sys.stderr)


if __name__ == "__main__":
    main()


