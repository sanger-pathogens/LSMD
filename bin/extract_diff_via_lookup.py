#!/usr/bin/env python3
"""
Compute a set difference (X - Y) via k-mer-exact lookup-based trimming,
avoiding the corrupted `sbwt difference` dummy-repair code path (see
notes/sbwt_difference_silent_corruption_bug.md).

Replaces the deleted extract_D_via_lookup.py, which kept a unitig only if
EVERY k-mer was absent from the target -- silently undercounting the true
difference at every unitig boundary that straddled a membership edge. This
version trims each unitig at the exact k-mer boundary instead: for every
maximal run of consecutive k-mer positions absent from the target (bitvector
'0'), it emits exactly the substring spanning that run. `sbwt build` then
re-derives exactly those k-mers from the emitted segments -- no gain, no
loss, mathematically equal to true set difference (see the reasoning in
notes/sbwt_difference_silent_corruption_bug.md for why this is exact, not
another approximation).

Input contract:
  --unitigs      Query unitigs FASTA (the "X" side), one header + one
                  sequence line per record (as produced by `sbwt dump-unitigs`
                  or Themisto2's export -- not line-wrapped).
  --bitvectors   Output of `sbwt lookup --index Y.sbwt --query <unitigs>
                  --membership-only`, one ASCII bitstring per unitig, in the
                  same order as --unitigs. B[i] = '1' iff the i-th k-mer of
                  that unitig is present in Y.

Prints kept/dropped k-mer counts to stderr in a machine-parseable
`KMER_COUNTS kept=<n> dropped=<n> total=<n>` line, so the caller can verify
correctness arithmetically after building:
    count(built result's k-mers) must equal `kept`
    `kept + dropped` must equal the source index's total k-mer count
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--unitigs", required=True, type=Path,
                         help="Query unitigs FASTA (the X side). Must be "
                              "unwrapped: exactly 1 header + 1 sequence line "
                              "per record.")
    parser.add_argument("--bitvectors", required=True, type=Path,
                         help="Output of `sbwt lookup --membership-only "
                              "--query <unitigs>` against Y, one ASCII "
                              "bitstring per line, same order as --unitigs.")
    parser.add_argument("--k", required=True, type=int,
                         help="k-mer length (must match the index build).")
    parser.add_argument("--out", required=True, type=Path,
                         help="Output FASTA path for the trimmed segments "
                              "(feed this to `sbwt build` next).")
    return parser.parse_args()


def main():
    args = parse_args()
    k = args.k

    n_unitigs = 0
    n_kept = 0      # k-mers absent from Y -> belong to X - Y
    n_dropped = 0   # k-mers present in Y -> excluded
    n_segments = 0

    with open(args.unitigs) as fa, open(args.bitvectors) as bv, open(args.out, "w") as out:
        while True:
            header = fa.readline()
            if not header:
                break
            seq = fa.readline().rstrip("\n")
            bits_line = bv.readline()

            if not header.startswith(">"):
                sys.exit(f"Expected FASTA header at unitig {n_unitigs}, got: {header!r}")
            if not seq:
                sys.exit(f"Unexpected EOF reading sequence for unitig {n_unitigs}")
            if not bits_line:
                sys.exit(f"Bitvector file ran out of lines at unitig {n_unitigs} "
                          f"(unitigs and bitvectors are out of sync)")

            bits = bits_line.rstrip("\n")
            expected_n_kmers = len(seq) - k + 1
            if expected_n_kmers < 1:
                sys.exit(f"Unitig {n_unitigs}: sequence length {len(seq)} is "
                          f"shorter than k={k}")
            if len(bits) != expected_n_kmers:
                sys.exit(f"Unitig {n_unitigs}: expected {expected_n_kmers} k-mer "
                          f"bit(s) for a sequence of length {len(seq)} (k={k}), "
                          f"got {len(bits)} -- unitigs/bitvectors out of sync "
                          f"or wrong k")

            # Walk the bitstring, emitting one FASTA record per maximal run
            # of consecutive '0' positions (k-mers absent from Y).
            run_start = None
            for i, ch in enumerate(bits):
                if ch == "1":
                    n_dropped += 1
                    if run_start is not None:
                        out.write(f">seg unitig={n_unitigs} start={run_start} end={i - 1}\n")
                        out.write(seq[run_start:i - 1 + k] + "\n")
                        n_segments += 1
                        run_start = None
                elif ch == "0":
                    n_kept += 1
                    if run_start is None:
                        run_start = i
                else:
                    sys.exit(f"Unitig {n_unitigs}, position {i}: unexpected "
                              f"bitvector character {ch!r} (expected '0' or '1')")

            # Flush a run that extends to the end of the unitig.
            if run_start is not None:
                last = len(bits) - 1
                out.write(f">seg unitig={n_unitigs} start={run_start} end={last}\n")
                out.write(seq[run_start:last + k] + "\n")
                n_segments += 1

            n_unitigs += 1

        # bitvector file must not have extra lines beyond the unitigs
        leftover = bv.readline()
        if leftover:
            sys.exit(f"Bitvector file has more lines than --unitigs has records "
                      f"(expected {n_unitigs}) -- unitigs and bitvectors are out of sync")

    print(f"Unitigs scanned    : {n_unitigs:,}", file=sys.stderr)
    print(f"Segments emitted   : {n_segments:,}", file=sys.stderr)
    print(f"K-mers kept        : {n_kept:,}", file=sys.stderr)
    print(f"K-mers dropped     : {n_dropped:,}", file=sys.stderr)
    print(f"KMER_COUNTS kept={n_kept} dropped={n_dropped} total={n_kept + n_dropped}", file=sys.stderr)


if __name__ == "__main__":
    main()
