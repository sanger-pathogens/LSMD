#!/usr/bin/env python3
"""
GENERALISED LINEAGE-CORE FILTERING PLAN
(species-agnostic, lineage-agnostic, background-db-agnostic)

Parameters (passed in, never hardcoded):
  SPECIES_INDEX     = combined SBWT/Themisto index for target species
                       (e.g. S. pneumoniae combined index, unitigs + colors)
  LINEAGE_LIST       = list of lineage IDs to process (e.g. all GPSCs of interest)
  LINEAGE_MAP        = mapping of lineage ID -> list of genome/color IDs
  BACKGROUND_DBS      = list of background reference indexes to check against
                       (e.g. [ATB, GTDB] — extensible to any future db)
  CORE_THRESHOLD      = core presence cutoff (e.g. 0.9 for "core", 0.0-0.1 for "catch-all")
  SPECIFICITY_CUTOFF  = max allowed presence outside lineage (e.g. 0.1)


STEP 0 — One-time setup (runs once, not per lineage)
  0.1  Load SPECIES_INDEX metadata (total genome/color count, color ID list)
  0.2  For each db in BACKGROUND_DBS:
         confirm db SBWT index exists / is built
         (if not, block — this is the GTDB-style dependency)
  0.3  Build A = species-wide k-mer set (already built once, reused for all lineages)
  0.4  For each db in BACKGROUND_DBS:
         F[db] = sbwt difference db.sbwt A.sbwt -o F_db.sbwt
         dump F[db] k-mers via sbwt dump-kmers  →  F_kmers[db]
         (these are species-agnostic backgrounds, built once, reused for every lineage)


STEP 1 — Per-lineage loop (repeat for each lineage in LINEAGE_LIST)

  for lineage in LINEAGE_LIST:

    1.1  Build lineage_membership boolean vector
         (implemented: build_colourid_to_lineage, find_lineage_colour_ids,
         build_lineage_membership_vector)
         lineage_membership = numpy.zeros(n_total_colors, dtype=bool)
         lineage_membership[LINEAGE_MAP[lineage]] = True

    1.2  (Optional, per lineage) Build B as real SBWT index
         — only if doing an SBWT-native D for benchmarking/validation
         — otherwise skip; D computed in Python at step 1.5

    1.3  Iterate SPECIES_INDEX unitigs one at a time
         (implemented: build_colorset_to_colourids, iter_unitig_lineage_stats)
         for each unitig:
           parse header → get color_set_id → look up colour set →
           build boolean vector `a`

           core_pct    = sum(a & lineage_membership) / sum(lineage_membership)
           outside_pct = sum(a & numpy.invert(lineage_membership)) / sum(numpy.invert(lineage_membership))

    1.4  Apply threshold → build C[lineage]
         (implemented: write_threshold_filtered_fasta)
         keep unitig if core_pct >= CORE_THRESHOLD
         (run twice per lineage if both core and catch-all modes are needed:
          CORE_THRESHOLD = 0.9-1.0 for core mode
          CORE_THRESHOLD = 0.0-0.1 for catch-all mode)
         write passing unitigs to C[lineage].fasta

    1.5  Compute D[lineage] = species-wide minus this lineage
         — SBWT-native (standard path):
             sbwt difference A.sbwt B_lineage.sbwt -o D_lineage.sbwt
             (queried directly via sbwt lookup in step 1.6, no dump-kmers needed)
         — Python-native version is BENCHMARK-ONLY (see GPSC1 benchmark ticket),
           only meaningful when compared at the threshold=100% boundary:
             for each unitig in SPECIES_INDEX:
               in_D = a & numpy.invert(lineage_membership)

    1.6  Compute E[lineage] = C[lineage] − D[lineage]     [NOT YET IMPLEMENTED]
         sbwt lookup -i D_lineage.sbwt -q C[lineage].fasta -o membership.txt --membership-only
         for each unitig in C[lineage], read its membership bitvector line:
           # TODO (decide when writing this step): exact keep/drop rule from
           # the bitvector — e.g. drop if ANY bit is 1 (strict), or drop if
           # >X% of bits are 1 (lenient). Not yet decided.
         write survivors → E[lineage].fasta

    1.7  Compute G[lineage][db] for each db in BACKGROUND_DBS   [NOT YET IMPLEMENTED]
         for db in BACKGROUND_DBS:
           sbwt lookup -i F_db.sbwt -q E[lineage].fasta -o membership.txt --membership-only
           for each unitig in E[lineage], read its membership bitvector line:
             # TODO: same bitvector rule decision as step 1.6, applied against F_db
           write survivors → G[lineage][db].fasta

    1.8  (Optional) apply SPECIFICITY_CUTOFF as an additional filter
         (implemented as part of write_threshold_filtered_fasta's
         specificity_cutoff parameter — applied alongside core_pct at step 1.4;
         re-applying as a final pass on G[lineage][db] is not yet implemented)
         keep only unitigs where outside_pct < SPECIFICITY_CUTOFF

    1.9  Write per-lineage summary log     [PARTIALLY IMPLEMENTED]
         (write_threshold_filtered_fasta writes total/kept/dropped counts for
         the C-building stage via stats_out_path; equivalent logging for
         E/G stages not yet implemented)


STEP 2 — Aggregate outputs
  2.1  Collate G[lineage][db] for all lineages × all backgrounds into final
       candidate marker set per lineage
  2.2  Reassemble into unitig FASTAs for LexicMap validation (existing PAT-3432)
  2.3  Compile summary stats across all lineages for sprint reporting


NOTES ON GENERALISATION
  - Nothing above references "S. pneumoniae," "GPSC," "ATB," or "GTDB" by name —
    all species/lineage/db specifics live in SPECIES_INDEX, LINEAGE_LIST,
    LINEAGE_MAP, and BACKGROUND_DBS, passed as parameters/config.
  - Adding a new background db (e.g. a third reference beyond ATB/GTDB) only
    requires adding it to BACKGROUND_DBS — steps 0.4 and 1.7 loop over it
    automatically, no code changes needed.
  - Adding a new target species requires only a new SPECIES_INDEX and its
    own LINEAGE_LIST/LINEAGE_MAP — same script, same logic.
  - Step 1.2 (per-lineage SBWT index) stays optional/config-driven, so the
    GPSC1 benchmark (PAT-3xxx) can validate Fork 1 vs Fork 2 without forcing
    every other lineage to pay the index-build cost.
  - `sbwt difference` does the heavy, index-native subtraction (F, F_gtdb, D)
    — this script does not reimplement that. For E and G, this script queries
    the resulting indexes directly via `sbwt lookup --membership-only` rather
    than loading a dumped k-mer list into memory — this matters at scale,
    since F (e.g. ATB minus species) could contain hundreds of millions of
    k-mers, and loading that into a Python set would reintroduce the same
    kind of memory blowup the original matrix-script redesign was meant to
    eliminate.

NOTE ON FILE SCOPE
  This file (lineage_core_filter.py) contains only the new lineage-core
  filtering logic. The existing species-exclusion script (used for PAT-3461,
  producing F via colour-based filtering) stays separate, unchanged, in
  filter_species_from_atb.py.
"""

import argparse
import gzip
import sys
from pathlib import Path
import numpy


def _open(path):
    """Open a file for reading, auto-detecting gzip regardless of the given path suffix.
    (Verbatim from build_kmer_matrix.py.)"""
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-lineage core presence and outside-lineage "
                     "specificity for unitigs in a Themisto2-exported species "
                     "index, then apply a threshold to build a filtered "
                     "candidate marker set (C).",
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
                         help="Output path for filtered unitig FASTA (C). "
                              "Use .gz suffix to compress.")
    parser.add_argument("--stats-out", type=Path, default=None,
                         help="Optional path to write summary stats "
                              "(total/kept/dropped counts).")
    return parser.parse_args()

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
    Stream export.unitigs.fa, applying threshold filtering to build C:
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
            f"Kept (C)              : {kept_count:,}",
            f"Dropped               : {dropped_count:,}",
            f"Core threshold        : {core_threshold}",
            f"Specificity cutoff    : {specificity_cutoff}",
        ]
        Path(stats_out_path).write_text("\n".join(stats_lines) + "\n")

    return total, kept_count, dropped_count


# TODO: E/G computation (steps 1.6/1.7) — blocked on bitvector collapse rule

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

    print(f"Done. Total: {total:,}  Kept (C): {kept_count:,}  "
          f"Dropped: {dropped_count:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
def parse_args():
    parser.add_argument(
        "-m", "--mode",
        choices=["core", "catchall"],
        default="core",
        help="'core': strict, present in 95%% of lineage genomes (override with --threshold). "
             "'catchall': relaxed, present in >50%% of lineage genomes (override with --threshold).",
    )
    parser.add_argument("-t", "--threshold", type=float, default=None,
        help="Override default cutoff (0.95 core / 0.5 catchall).")
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
