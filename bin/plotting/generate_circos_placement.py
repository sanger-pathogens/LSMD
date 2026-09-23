#!/usr/bin/env python3
"""
Marker-placement Circos figure: one anchor genome (the target lineage's own
reference, e.g. c6706 for 7PET), every candidate unitig placed on it by blastn,
tiled and coloured PASS (cleared the length/GC filter) vs REJECT (too short/
off-GC). Within-genome links then connect markers sitting within
--merge-distance bp of each other on the same contig -- candidate clusters
where a primer pair could span two-or-more short REJECT unitigs plus their
(non-specific) intervening sequence as one longer amplicon, rather than
needing any single unitig to individually clear the length floor.

This answers a different question from circos_tier3 (generate_circos_tier3.py):
that figure links a marker's position on the anchor to its position in a
DIFFERENT lineage's genome (specificity: "does this also occur elsewhere?").
This one only ever looks at the anchor genome -- it's asking "where do our own
candidates sit, and which short ones are close enough to merge?"

Runs its own blastn (markers vs the anchor only -- unlike the population run
tier3 reuses, no pre-existing hit table covers this).

Usage:
    python3 generate_circos_placement.py <pass.fasta> <reject.fasta> <anchor.fasta> \\
        [--merge-distance 300] [--anchor-label c6706] [-o circos_placement] \\
        [--etc-dir circos_tier3/etc] [--render]
"""
import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PASS_COLOR = "blue"  # matches marker_heatmap.py / group_distribution.py: blue = the on-target story
REJECT_COLOR = (
    "grey"  # grey = not (yet) usable, not "wrong" -- same convention as group_distribution.py's "everything else"
)
LINK_COLOR = "dgreen"  # candidate-merge links: a distinct third colour, its own channel (status), not magnitude

BLAST_FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovhsp evalue bitscore"

IDEOGRAM_CONF = """<ideogram>
<spacing>
default = 0.005r
</spacing>
radius    = 0.85r
thickness = 60p
fill      = yes
show_label = yes
label_font = default
label_radius = dims(ideogram,radius_outer) + 30p
label_size = 28
label_parallel = yes
</ideogram>
"""

TICKS_CONF = """show_ticks       = yes
show_tick_labels = yes
<ticks>
radius      = 1r
color       = black
thickness   = 2p
multiplier  = 1e-6
format      = %d

<tick>
spacing     = 0.1u
size        = 5p
</tick>
<tick>
spacing     = 0.5u
size        = 10p
show_label  = yes
label_size  = 16p
label_offset = 5p
suffix      = " Mb"
</tick>
</ticks>
"""


def sanitize(genome_id):
    return re.sub(r"[^A-Za-z0-9_]", "_", genome_id)


def run_blastn(query_fasta, anchor_fasta, outdir):
    db_dir = Path(outdir) / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "anchor"
    # no -parse_seqids: it wraps sseqid as "ref|ACCESSION|" in blastn's own
    # output, which then wouldn't match the bare accession .fai/karyotype use.
    r = subprocess.run(
        ["makeblastdb", "-in", str(anchor_fasta), "-dbtype", "nucl", "-out", str(db_prefix)], capture_output=True
    )
    if r.returncode != 0:
        sys.exit(f"makeblastdb failed:\nstdout: {r.stdout.decode()}\nstderr: {r.stderr.decode()}")

    hits_path = Path(outdir) / "markers_vs_anchor.tsv"
    with open(hits_path, "w") as out:
        subprocess.run(
            [
                "blastn",
                "-query",
                str(query_fasta),
                "-db",
                str(db_prefix),
                "-task",
                "blastn",
                "-word_size",
                "11",
                "-dust",
                "no",
                "-perc_identity",
                "90",
                "-evalue",
                "1",
                "-max_target_seqs",
                "5",
                "-outfmt",
                BLAST_FMT,
            ],
            check=True,
            stdout=out,
        )
    return hits_path


def load_hits(hits_path):
    """qseqid -> best hit (highest bitscore) as dict."""
    best = {}
    with open(hits_path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, qcovhsp, evalue, bitscore = (
                cols
            )
            bitscore = float(bitscore)
            if qseqid not in best or bitscore > best[qseqid]["bitscore"]:
                sstart, send = int(sstart), int(send)
                # -parse_seqids wraps sseqid as e.g. "ref|NZ_CP064351.1|" -- .fai
                # (and so our karyotype contig keys) use the bare accession.
                contig = sseqid.split("|")[1] if sseqid.count("|") >= 2 else sseqid
                best[qseqid] = {
                    "contig": contig,
                    "start": min(sstart, send),
                    "end": max(sstart, send),
                    "bitscore": bitscore,
                    "pident": float(pident),
                }
    return best


def contig_lengths(fasta_path, index_dir):
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    fai = index_dir / (Path(fasta_path).name + ".fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(fasta_path), "--fai-idx", str(fai)], check=True)
    lengths = {}
    with open(fai) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def build_karyotype(anchor_label, lengths, used_contigs, outdir):
    lines = []
    contig_to_chr = {}
    contigs = [c for c in lengths if c in used_contigs]
    multi = len(contigs) > 1
    for i, contig in enumerate(contigs):
        chr_id = f"{sanitize(anchor_label)}_{i}"
        contig_to_chr[contig] = chr_id
        # distinguish replicons in the label when there's more than one --
        # otherwise every ideogram segment prints the identical "c6706" label
        # on top of itself. Karyotype labels can't contain whitespace (circos
        # splits the line on it), hence "_" rather than " (i/n)".
        label = f"{anchor_label}_{i + 1}of{len(contigs)}" if multi else anchor_label
        lines.append(f"chr - {chr_id} {label} 0 {lengths[contig]} vdgrey")
    (Path(outdir) / "karyotype.txt").write_text("\n".join(lines) + "\n")
    print(f"wrote {outdir}/karyotype.txt ({len(lines)} replicons)", file=sys.stderr)
    return contig_to_chr


# Markers are 31-1311bp on a ~4Mb genome -- at true scale most are sub-pixel.
# Padded by this many bp on each side purely so the tile is visible; it's a
# rendering aid, not a claim about the marker's real footprint (the real
# start/end feed links.txt and the merge-distance calculation untouched).
VISUAL_PAD_BP = 4000


def build_markers_track(placements, contig_to_chr, outdir):
    """One line per placed marker: <chr> <start> <end> id=<marker_id> fill_color=<...>"""
    lines = []
    for marker_id, verdict, hit in placements:
        chr_id = contig_to_chr.get(hit["contig"])
        if chr_id is None:
            continue
        color = PASS_COLOR if verdict == "PASS" else REJECT_COLOR
        start = max(0, hit["start"] - VISUAL_PAD_BP)
        end = hit["end"] + VISUAL_PAD_BP
        lines.append(f"{chr_id} {start} {end} id={marker_id} fill_color={color}")
    (Path(outdir) / "markers.txt").write_text("\n".join(lines) + "\n")
    print(f"wrote {outdir}/markers.txt ({len(lines)} placed markers)", file=sys.stderr)


def build_merge_links(placements, contig_to_chr, merge_distance, outdir):
    """Any two markers on the same contig within merge_distance bp of each
    other -> one link. Only REJECT markers matter for the "could this become
    usable" story, but a REJECT next to a PASS is worth showing too (it may
    let you extend an already-working amplicon)."""
    by_contig = defaultdict(list)
    for marker_id, verdict, hit in placements:
        chr_id = contig_to_chr.get(hit["contig"])
        if chr_id is None:
            continue
        by_contig[chr_id].append((hit["start"], hit["end"], marker_id, verdict))

    lines = []
    n_links = 0
    n_clustered_reject = set()
    for chr_id, entries in by_contig.items():
        entries.sort()
        for i in range(len(entries)):
            s1, e1, id1, v1 = entries[i]
            for j in range(i + 1, len(entries)):
                s2, e2, id2, v2 = entries[j]
                gap = s2 - e1
                if gap > merge_distance:
                    break  # sorted by start -- nothing further on this contig is closer
                lines.append(f"{chr_id} {s1} {e1} {chr_id} {s2} {e2} color={LINK_COLOR}")
                n_links += 1
                if v1 == "REJECT":
                    n_clustered_reject.add(id1)
                if v2 == "REJECT":
                    n_clustered_reject.add(id2)

    (Path(outdir) / "links.txt").write_text("\n".join(lines) + "\n")
    print(
        f"wrote {outdir}/links.txt ({n_links} candidate-merge links within {merge_distance}bp; "
        f"{len(n_clustered_reject)} otherwise-rejected markers sit in a cluster)",
        file=sys.stderr,
    )
    return n_clustered_reject


def build_conf(outdir, etc_dir, image_name):
    etc_rel = Path(etc_dir).name if Path(etc_dir).parent == Path(outdir) else str(etc_dir)
    conf = f"""# Generated by generate_circos_placement.py.
karyotype         = karyotype.txt
chromosomes_units = 1000000

<<include ideogram.conf>>
<<include ticks.conf>>

<highlights>
<highlight>
file       = markers.txt
r0         = 0.82r
r1         = 0.99r
</highlight>
</highlights>

<links>
<link>
file          = links.txt
radius        = 0.85r
bezier_radius = 0.1r
thickness     = 2
color         = {LINK_COLOR}
</link>
</links>

<image>
dir               = .
file              = {image_name}
png               = yes
svg               = yes
radius            = 900p
background        = white
angle_offset      = -90
auto_alpha_colors = yes
auto_alpha_steps  = 5
</image>

<colors>
<<include {etc_rel}/colors.conf>>
</colors>
<fonts>
<<include {etc_rel}/fonts.conf>>
</fonts>
<patterns>
<<include {etc_rel}/patterns.conf>>
</patterns>
<<include {etc_rel}/housekeeping.conf>>
data_out_of_range* = trim
"""
    (Path(outdir) / "circos.conf").write_text(conf)
    (Path(outdir) / "ideogram.conf").write_text(IDEOGRAM_CONF)
    (Path(outdir) / "ticks.conf").write_text(TICKS_CONF)
    print(
        f"wrote {outdir}/circos.conf (+ ideogram.conf, ticks.conf) -- blue tile = PASS, grey tile = REJECT, "
        f"{LINK_COLOR} link = within merge distance",
        file=sys.stderr,
    )


def load_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pass_fasta", help="post-processing PASS markers (filtered_unitigs.fasta)")
    p.add_argument("reject_fasta", help="post-processing REJECT markers (rejected_unitigs.fasta)")
    p.add_argument("anchor_fasta", help="one reference genome the target lineage carries (e.g. c6706)")
    p.add_argument("--anchor-label", default="c6706")
    p.add_argument(
        "--merge-distance",
        type=int,
        default=300,
        help="max bp gap between two markers to draw a candidate-merge link (default: 300)",
    )
    p.add_argument("-o", "--outdir", default="circos_placement", help="Output directory (default: circos_placement)")
    p.add_argument("--etc-dir", default=None, help="Existing circos etc/ dir (e.g. circos_tier3/etc)")
    p.add_argument("--image-name", default="circos_placement")
    p.add_argument("--render", action="store_true")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    query_fasta = outdir / "markers_query.fasta"
    with open(query_fasta, "w") as out:
        for path in (args.pass_fasta, args.reject_fasta):
            out.write(Path(path).read_text())

    hits_path = run_blastn(query_fasta, args.anchor_fasta, outdir)
    hits = load_hits(hits_path)

    pass_ids = set(load_fasta_ids(args.pass_fasta))
    reject_ids = set(load_fasta_ids(args.reject_fasta))
    placements = []
    n_unplaced = 0
    for marker_id, hit in hits.items():
        verdict = "PASS" if marker_id in pass_ids else ("REJECT" if marker_id in reject_ids else None)
        if verdict is None:
            continue
        placements.append((marker_id, verdict, hit))
    n_unplaced = len(pass_ids | reject_ids) - len(placements)
    print(
        f"placed {len(placements)} / {len(pass_ids) + len(reject_ids)} markers on {args.anchor_label} "
        f"({n_unplaced} had no blastn hit on the anchor)",
        file=sys.stderr,
    )

    lengths = contig_lengths(args.anchor_fasta, outdir / ".fai_cache")
    used_contigs = {hit["contig"] for _, _, hit in placements}
    contig_to_chr = build_karyotype(args.anchor_label, lengths, used_contigs, outdir)
    build_markers_track(placements, contig_to_chr, outdir)
    clustered_reject = build_merge_links(placements, contig_to_chr, args.merge_distance, outdir)

    etc_dir = args.etc_dir
    if etc_dir is None or not Path(etc_dir).exists():
        sys.exit(f"--etc-dir not found ({etc_dir}) -- pass an existing circos etc/ (e.g. circos_tier3/etc)")
    build_conf(outdir, etc_dir, args.image_name)

    print(
        f"\n{len(clustered_reject)} of {len(reject_ids)} REJECTed (too-short) markers sit within "
        f"{args.merge_distance}bp of another marker -- worth checking those clusters for a mergeable amplicon.",
        file=sys.stderr,
    )

    if args.render:
        subprocess.run(["circos", "-conf", "circos.conf"], cwd=outdir, check=True)
        print(f"rendered {outdir}/{args.image_name}.png / .svg", file=sys.stderr)
    else:
        print(f"run: cd {outdir} && circos -conf circos.conf", file=sys.stderr)


if __name__ == "__main__":
    main()
