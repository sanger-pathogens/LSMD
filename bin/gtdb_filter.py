#!/usr/bin/env python3
"""
GTDB specificity filter for candidate marker unitigs.

For each unitig, a PASS requires:
  - zero GTDB hits, OR
  - hits ONLY to genomes classified as s__Streptococcus pneumoniae in GTDB

A FAIL occurs if any hit maps to a non-S. pneumoniae genome.

Steps:
  1. Build a set of S. pneumoniae accessions from bac120_taxonomy_r226.tsv
  2. Run minimap2 of unitigs against the GTDB genome collection
     (done externally via the companion bsub script)
  3. Parse PAF output and classify each unitig
  4. Write passing and failing FASTAs

Usage:
  python gtdb_filter.py --candidates gpsc1_strict.fa --paf hits.paf --output gtdb_pass.fa
"""

import argparse
import os
import sys
from collections import defaultdict

TAXONOMY = (
    "/data/pam/collections/GTDB/release226/genomic_files_all_minus_suppressed"
    "/metadata/2025-03-05/bac120_taxonomy_r226.tsv"
)

PNEUMO_SPECIES = "s__Streptococcus pneumoniae"


def load_spneumo_accessions(taxonomy_path):
    """Return set of accession IDs (stripped of RS_/GB_ prefix) that are S. pneumoniae."""
    spneumo = set()
    with open(taxonomy_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            acc, lineage = line.split('\t', 1)
            # strip RS_ or GB_ prefix to get bare accession
            bare = acc[3:] if acc.startswith(('RS_', 'GB_')) else acc
            if PNEUMO_SPECIES in lineage:
                spneumo.add(bare)
    return spneumo


def parse_paf(paf_path, spneumo_accessions):
    """
    Parse minimap2 PAF output.
    Target name in PAF is the genome filename (stem), e.g. GCA_000435495.1
    Returns:
      pass_ids : set of query names with zero non-pneumo hits
      fail_ids : set of query names with >=1 non-pneumo hit
    """
    hit_status = defaultdict(lambda: "no_hit")  # query -> "pass" | "fail"

    with open(paf_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            query = parts[0]
            target = parts[5]  # target sequence name

            # target name from minimap2 is the fasta header; we expect accession as prefix
            # e.g. GCA_000435495.1 or the filename stem
            acc = target.split()[0]
            # strip path if accidentally included
            acc = os.path.basename(acc)
            # strip _genomic suffix and .fna
            acc = acc.replace('_genomic.fna', '').replace('.fna', '')

            if acc in spneumo_accessions:
                if hit_status[query] == "no_hit":
                    hit_status[query] = "pass"
                # already pass or fail — pneumo hit doesn't change fail
            else:
                hit_status[query] = "fail"

    pass_ids = {q for q, s in hit_status.items() if s in ("pass", "no_hit")}
    fail_ids = {q for q, s in hit_status.items() if s == "fail"}
    return pass_ids, fail_ids


def filter_fasta(input_fa, pass_ids, pass_out, fail_out):
    """Split input FASTA into pass and fail based on header name."""
    n_pass = n_fail = n_no_hit = 0
    pass_f = open(pass_out, 'w')
    fail_f = open(fail_out, 'w')
    try:
        with open(input_fa) as f:
            header = None
            seq = None
            name = None
            for line in f:
                line = line.rstrip()
                if line.startswith('>'):
                    if header and seq:
                        if name in fail_ids:
                            fail_f.write(f"{header}\n{seq}\n")
                            n_fail += 1
                        else:
                            pass_f.write(f"{header}\n{seq}\n")
                            if name in pass_ids:
                                n_pass += 1
                            else:
                                n_no_hit += 1
                    header = line
                    name = line[1:].split()[0]
                    seq = None
                else:
                    seq = line
            if header and seq:
                if name in fail_ids:
                    fail_f.write(f"{header}\n{seq}\n")
                    n_fail += 1
                else:
                    pass_f.write(f"{header}\n{seq}\n")
                    if name in pass_ids:
                        n_pass += 1
                    else:
                        n_no_hit += 1
    finally:
        pass_f.close()
        fail_f.close()

    return n_pass, n_no_hit, n_fail


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidates", "-c", required=True,
                        help="Input FASTA of strict-marker unitigs (from gpsc_stats_extract.py)")
    parser.add_argument("--paf", required=True,
                        help="PAF file from minimap2 search against GTDB")
    parser.add_argument("--output", "-o", required=True,
                        help="Output prefix; writes <prefix>_gtdb_pass.fa and <prefix>_gtdb_fail.fa")
    parser.add_argument("--taxonomy", default=TAXONOMY,
                        help="GTDB bac120 taxonomy TSV")
    args = parser.parse_args()

    print(f"Loading S. pneumoniae accessions from taxonomy...", file=sys.stderr)
    spneumo = load_spneumo_accessions(args.taxonomy)
    print(f"  {len(spneumo)} S. pneumoniae accessions", file=sys.stderr)

    print(f"Parsing PAF hits...", file=sys.stderr)
    pass_ids, fail_ids = parse_paf(args.paf, spneumo)
    print(f"  {len(pass_ids)} queries with 0 or only S. pneumoniae hits", file=sys.stderr)
    print(f"  {len(fail_ids)} queries with non-S. pneumoniae hits", file=sys.stderr)

    pass_out = f"{args.output}_gtdb_pass.fa"
    fail_out = f"{args.output}_gtdb_fail.fa"
    n_pass, n_no_hit, n_fail = filter_fasta(args.candidates, pass_ids, fail_ids, pass_out, fail_out)

    print(f"\nResults:", file=sys.stderr)
    print(f"  No GTDB hit:          {n_no_hit}", file=sys.stderr)
    print(f"  S. pneumoniae hit only: {n_pass}", file=sys.stderr)
    print(f"  Non-S. pneumoniae hit:  {n_fail} (written to {fail_out})", file=sys.stderr)
    print(f"  PASS total: {n_no_hit + n_pass} → {pass_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
