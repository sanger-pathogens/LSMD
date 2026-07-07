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

    1.1  Build lineagecore boolean vector
         lineagecore = numpy.zeros(n_total_colors, dtype=bool)
         lineagecore[LINEAGE_MAP[lineage]] = True

    1.2  (Optional, per lineage) Build B as real SBWT index
         — only if doing an SBWT-native D for benchmarking/validation
         — otherwise skip; D computed in Python at step 1.5

    1.3  Iterate SPECIES_INDEX unitigs one at a time
         for each unitig:
           parse header → get color set → build boolean vector `a`

           core_pct    = sum(a & lineagecore) / sum(lineagecore)
           outside_pct = sum(a & numpy.invert(lineagecore)) / sum(numpy.invert(lineagecore))

    1.4  Apply threshold → build C[lineage]
         keep unitig if core_pct >= CORE_THRESHOLD
         (run twice per lineage if both core and catch-all modes are needed:
          CORE_THRESHOLD = 0.9-1.0 for core mode
          CORE_THRESHOLD = 0.0-0.1 for catch-all mode)
         write passing unitigs to C[lineage].fasta

    1.5  Compute D[lineage] = species-wide minus this lineage
         — Python-native version (no per-lineage index needed):
             for each unitig in SPECIES_INDEX:
               in_D = a & numpy.invert(lineagecore)
         — OR SBWT-native version if B[lineage] was built (step 1.2):
             sbwt difference A.sbwt B_lineage.sbwt -o D_lineage.sbwt
             dump D_lineage k-mers → D_kmers[lineage]

    1.6  Compute E[lineage] = C[lineage] − D[lineage]
         for each unitig in C[lineage]:
           drop if its k-mers appear in D_kmers[lineage] (or in-memory D from 1.5)
         write survivors → E[lineage].fasta

    1.7  Compute G[lineage][db] for each db in BACKGROUND_DBS
         for db in BACKGROUND_DBS:
           for each unitig in E[lineage]:
             drop if its k-mers appear in F_kmers[db]
           write survivors → G[lineage][db].fasta

    1.8  (Optional) apply SPECIFICITY_CUTOFF as an additional filter
         keep only unitigs where outside_pct < SPECIFICITY_CUTOFF
         — can be applied at step 1.4 alongside core_pct, or as a final
           pass on G[lineage][db]

    1.9  Write per-lineage summary log
         (unitig counts at each stage: species-wide → C → E → G[db] per db,
          plus core_pct/outside_pct distribution for QC)


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
"""

import argparse
import gzip
import re
import sys
from pathlib import Path


# numpy boolean vectors:
"""
What &, invert(), sum() do
and explain whyt each is used in their areas etc
"""

# sbwt lookup --membership-only
"""
explore what this comand actually returns.
A bitvector (one bit per k-mer in your query sequence) in k-mer order along the unitig and not per unitig, per k-mer
    decide how to collapse that bitvactor into a single keep/drop decision per unitig
    understand this before passig the loop as it changes what you do with th eoutput
"""

# where color IDS come from and what they map to
"""
Colour IDs in Themisto's colour space aren't inherently meaningful — they only become "genome X" or "lineage Y" via color_names.txt. You need to be clear on the chain: colour_id → filepath → genome accession → GPSC assignment (your existing lookup table). This is currently species-level (find_species_colour_ids) and needs converting to lineage-level — make sure you know exactly where your GPSC-to-genome mapping lives and what format it's in before you start writing find_lineage_colour_ids().
"""


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
        description="Exclude a species' k-mers/unitigs from a Themisto-exported background set.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--unitigs", required=True, type=Path,
                         help="export.unitigs.fa from Themisto export.")
    parser.add_argument("--color-sets", required=True, type=Path,
                         help="export.color_sets.txt(.gz) from Themisto export.")
    parser.add_argument("--colour-map", required=True, type=Path,
                         help="color_names.txt — colour ID to species/genome name mapping.")
    parser.add_argument("--species", required=True,
                         help='Species to exclude, e.g. "Streptococcus pneumoniae". '
                              "Matched as a case-insensitive substring against colour_map entries.")
    parser.add_argument("--out", required=True, type=Path,
                         help="Output path for filtered unitig FASTA (F). Use .gz suffix to compress.")
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

def find_species_colour_ids(colour_map_path, species):
    """
    Parse color_names.txt and return the set of colour IDs matching `species`.

    color_names.txt lines: "<colour_id>\\t<path>", path like
    "per_species_unitigs/<genus>_<species>[suffix]-unitigs-k31.fna". Suffix is
    either an accession/strain tag ("_sp011388195") or a bare shard letter
    ("streptococcus_pneumoniaea") when ATB splits an oversampled species
    across multiple colour files.

    Normalises the species name to "genus_species" and matches it against the
    filename stem, allowing either suffix form.
    """
    # normalise "Streptococcus pneumoniae" -> "streptococcus_pneumoniae"
    normalised_species = "_".join(species.lower().split())
    # exact match, OR "<species>_<accession/strain>", OR "<species><shard letter>"
    stem_pattern = re.compile(
        rf"^{re.escape(normalised_species)}([a-z]|_.+)?$"
    )
    matched_ids = set()

    with open(colour_map_path) as fh:
        for line_num, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split(None, 1)
            if len(parts) < 2:
                print(f"Warning: could not parse line {line_num} in colour map: {line!r}",
                      file=sys.stderr)
                continue
            colour_id_str, path = parts[0], parts[1]

            filename = path.rsplit("/", 1)[-1].lower()  # drop "per_species_unitigs/"
            # strip the trailing "-unitigs-k31.fna" (or similar) suffix
            stem = filename.split("-unitigs")[0]

            if stem_pattern.match(stem):
                matched_ids.add(int(colour_id_str))

    if not matched_ids:
        raise ValueError(
            f"No colour IDs matched species '{species}' in {colour_map_path}. "
            "Check the species name/spelling and the colour map file format."
        )

    print(f"Matched {len(matched_ids)} colour ID(s) for '{species}': "
          f"{sorted(matched_ids)[:10]}{'...' if len(matched_ids) > 10 else ''}",
          file=sys.stderr)
    return matched_ids


def find_excluded_colorsets(color_sets_path, species_colour_ids):
    """
    Stream export.color_sets.txt (format: "color_set_id=N size=M c1 c2 c3 ...").
    Returns the set of color_set_ids whose colour list intersects
    species_colour_ids — these are the colour sets to exclude.
    """
    excluded_colorsets = set()
    with _open(color_sets_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            color_set_id = int(parts[0].split("=")[1])
            colours = {int(tok) for tok in parts[2:]}  # skip color_set_id=N, size=M
            if colours & species_colour_ids:
                excluded_colorsets.add(color_set_id)
    return excluded_colorsets


# iterate unitigs one at a time mechanism needed. reuse this loop structre, change what happens inside flush()
# — instead of exclude/retain by colorset membership, you'll be building the boolean vector a and computing core_pct/outside_pct.
def write_filtered_fasta(unitigs_path, excluded_colorsets, out_path):
    """
    Single streaming pass over export.unitigs.fa, writing out only unitigs
    whose colour set is NOT in excluded_colorsets.

    Each unitig header already contains its own color_set_id
    (">unitig_id=N color_set_id=M"), so no separate lookup pass over
    unitigs.fa is needed before this one — the filter decision is made
    directly from the header as we stream through.

    Returns (total_unitigs, excluded_count, retained_count).
    """
    total = 0
    excluded_count = 0
    retained_count = 0

    with _open(unitigs_path) as in_fh, _open_out(out_path) as out_fh:
        current_header = None
        current_colorset = None
        current_seq = None

        def flush():
            nonlocal excluded_count, retained_count
            if current_header is None:
                return
            if current_colorset in excluded_colorsets:
                excluded_count += 1
            else:
                out_fh.write(f">{current_header}\n{current_seq}\n")
                retained_count += 1

        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                total += 1
                current_header = line[1:].strip()
                current_colorset = None
                for tok in current_header.split():
                    if tok.startswith("color_set_id="):
                        current_colorset = int(tok.split("=", 1)[1])
                        break
                current_seq = None
            else:
                current_seq = line
        flush()  # last record

    return total, excluded_count, retained_count


def main():
    args = parse_args()

    print(f"Looking up colour ID(s) for species: {args.species}", file=sys.stderr)
    species_colour_ids = find_species_colour_ids(args.colour_map, args.species)

    print("Scanning colour sets for species matches...", file=sys.stderr)
    excluded_colorsets = find_excluded_colorsets(args.color_sets, species_colour_ids)
    print(f"  {len(excluded_colorsets):,} colour sets contain '{args.species}'", file=sys.stderr)

    print(f"Writing filtered FASTA to {args.out} ...", file=sys.stderr)
    total, excluded_count, retained_count = write_filtered_fasta(
        args.unitigs, excluded_colorsets, args.out
    )

    stats_path = Path(str(args.out).split(".fa")[0] + ".stats.txt")
    stats_lines = [
        f"Species excluded        : {args.species}",
        f"Colour IDs matched       : {sorted(species_colour_ids)}",
        f"Total unitigs scanned    : {total:,}",
        f"Unitigs excluded         : {excluded_count:,}",
        f"Unitigs retained (F)     : {retained_count:,}",
    ]
    stats_path.write_text("\n".join(stats_lines) + "\n")
    for line in stats_lines:
        print(line)


if __name__ == "__main__":
    main()