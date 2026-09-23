#!/usr/bin/env python3
"""
Reproduce the (hand-built) circos_tier3 figure -- karyotype.txt + links.txt +
circos.conf/ideogram.conf/ticks.conf -- from a Tier-3 blastn hsp_detail.tsv,
so it's a repeatable build instead of a one-off hand edit.

What the figure shows: the anchor genome (a genome carrying the target group,
e.g. c6706 for 7PET) as the hub, plus one representative genome per *other*
group the markers also land in. A ribbon joins a marker's position on the
anchor to its position in that representative genome -- a dense fan of ribbons
to one lineage means the marker set is not target-specific w.r.t. that lineage.

Only one blastn run is needed, not two: hsp_detail.tsv already has the anchor's
own coordinates as just another row (genome == anchor), since the anchor is
itself one of the assemblies in the population blastn run. Anchor-vs-marker
coordinates and representative-vs-marker coordinates are joined on marker_id.

Representative-genome selection: the original circos_tier3/ was hand-built with
manually chosen reference strains (not reproducible from the table alone). This
script instead picks, per non-target group, the genome with the MOST marker
hits in that group -- a deterministic, scriptable stand-in for "a genome that
represents this lineage well" that tells the same story (dense fan = shared).
Pass --representatives to pin exact genome IDs instead (e.g. to match the
original figure's choices for a side-by-side comparison).

Usage:
    python3 generate_circos_tier3.py <hsp_detail.tsv> <fasta_dir> \\
        [--anchor c6706] [--target-group 7PET] [--n-lineages 8] \\
        [--representatives 10432-62=Non-7PET,mp_070116=Vibrioparacholerae,...] \\
        [--outdir circos_tier3] [--etc-dir circos_tier3/etc] [--render]
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Circos built-in colour names, fixed order -- matches the original figure's
# palette and follows the "assign categorical hues in fixed order, never
# cycled" rule. A 9th lineage folds into "grey" rather than generating a new hue.
LINEAGE_COLORS = ["blue", "purple", "orange", "green", "red", "vdblue", "dgreen", "brown"]
ANCHOR_COLOR = "vdgrey"
OTHER_COLOR = "grey"

IDEOGRAM_CONF = """<ideogram>
<spacing>
default = 0.005r
</spacing>
radius    = 0.80r
thickness = 40p
fill      = yes
show_label = yes
label_font = default
label_radius = dims(ideogram,radius_outer) + 30p
label_size = 24
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
spacing     = 0.5u
size        = 5p
</tick>
<tick>
spacing     = 1u
size        = 10p
show_label  = yes
label_size  = 16p
label_offset = 5p
suffix      = " Mb"
</tick>
</ticks>
"""


def load_hsp_detail(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def pick_representatives(rows, anchor, target_group, n_lineages, pinned):
    """group -> representative genome id."""
    if pinned:
        return dict(pinned)
    hits_by_group_genome = Counter()
    for r in rows:
        if r["genome"] == anchor or r["group"] in (target_group, "__unmapped__", ""):
            continue
        hits_by_group_genome[(r["group"], r["genome"])] += 1

    best_per_group = {}
    for (group, genome), n in hits_by_group_genome.items():
        if group not in best_per_group or n > best_per_group[group][1]:
            best_per_group[group] = (genome, n)

    ordered_groups = sorted(best_per_group, key=lambda g: -best_per_group[g][1])[:n_lineages]
    return {g: best_per_group[g][0] for g in ordered_groups}


def find_fasta(fasta_dir, genome_id):
    """Genome IDs in hsp_detail.tsv aren't guaranteed to match on-disk case
    (e.g. 'c6706' in the table vs 'C6706.fasta' on disk)."""
    d = Path(fasta_dir)
    for candidate in (f"{genome_id}.fasta", f"{genome_id}.fasta.gz"):
        p = d / candidate
        if p.exists():
            return p
    lower = genome_id.lower()
    for p in d.iterdir():
        if p.stem.lower() == lower or p.name.lower().startswith(lower + "."):
            return p
    return None


def contig_lengths(fasta_path, index_dir):
    """samtools faidx -> {contig: length}, in file order. The index is written
    under index_dir (not next to fasta_path) since source FASTA directories
    here are often read-only to this user (e.g. team216's vibriowatch dir)."""
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


def sanitize(genome_id):
    """Circos chromosome IDs can't contain '.' or '-' cleanly in this context
    (same convention the original karyotype.txt used, e.g. 2012env-9 -> 2012env_9)."""
    return genome_id.replace("-", "_").replace(".", "_")


def build_karyotype(genomes, fasta_dir, outdir, used_contigs):
    """genomes: [(genome_id, label, color)], anchor first. used_contigs:
    {genome_id: set(contig)} -- restricts the ideogram to contigs that
    actually carry a marker hit (draft assemblies here run to hundreds of
    contigs; the original hand-built figure only ever showed the handful with
    real data, keeping each genome to ~2-3 replicons instead of the whole draft).
    Returns {genome_id: {contig: chr_id}}."""
    lines = []
    contig_to_chr = {}
    for genome_id, label, color in genomes:
        fasta = find_fasta(fasta_dir, genome_id)
        if fasta is None:
            sys.exit(f"could not find a FASTA for genome '{genome_id}' under {fasta_dir}")
        lengths = contig_lengths(fasta, Path(outdir) / ".fai_cache")
        keep = used_contigs.get(genome_id, set(lengths))
        contig_to_chr[genome_id] = {}
        # karyotype labels can't contain whitespace (circos splits the line
        # on it) -- group names here include e.g. "Gulf Coast", "ELA-3 / part of L9".
        safe_label = re.sub(r"\s+", "_", label)
        for i, contig in enumerate(c for c in lengths if c in keep):
            length = lengths[contig]
            chr_id = f"{sanitize(genome_id)}_{i}"
            contig_to_chr[genome_id][contig] = chr_id
            lines.append(f"chr - {chr_id} {safe_label} 0 {length} {color}")

    out_path = Path(outdir) / "karyotype.txt"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_path} ({len(lines)} replicons)", file=sys.stderr)
    return contig_to_chr


def build_links(rows, anchor, representatives, contig_to_chr, outdir):
    """One ribbon per marker shared between the anchor and a representative
    genome, joined on marker_id."""
    anchor_hits = defaultdict(list)
    for r in rows:
        if r["genome"] == anchor:
            anchor_hits[r["marker_id"]].append(r)

    rep_genomes = set(representatives.values())
    rep_hits = defaultdict(list)
    for r in rows:
        if r["genome"] in rep_genomes:
            rep_hits[(r["genome"], r["marker_id"])].append(r)

    lines = []
    color_by_group = dict(zip(representatives.keys(), LINEAGE_COLORS))
    n_ribbons = 0
    for group, genome_id in representatives.items():
        color = color_by_group[group]
        for marker_id, a_rows in anchor_hits.items():
            g_rows = rep_hits.get((genome_id, marker_id))
            if not g_rows:
                continue
            for a in a_rows:
                a_chr = contig_to_chr[anchor].get(a["contig"])
                if a_chr is None:
                    continue
                a_start, a_end = sorted((int(a["sstart"]), int(a["send"])))
                for g in g_rows:
                    g_chr = contig_to_chr[genome_id].get(g["contig"])
                    if g_chr is None:
                        continue
                    g_start, g_end = sorted((int(g["sstart"]), int(g["send"])))
                    lines.append(f"{a_chr} {a_start} {a_end} {g_chr} {g_start} {g_end} color={color}")
                    n_ribbons += 1

    out_path = Path(outdir) / "links.txt"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_path} ({n_ribbons} ribbons across {len(representatives)} lineages)", file=sys.stderr)


def build_conf(outdir, etc_dir, image_name):
    etc_rel = Path(etc_dir).name if Path(etc_dir).parent == Path(outdir) else etc_dir
    conf = f"""# Generated by generate_circos_tier3.py -- see module docstring for how this
# reproduces the hand-built circos_tier3/ figure.
karyotype         = karyotype.txt
chromosomes_units = 1000000

<<include ideogram.conf>>
<<include ticks.conf>>

<links>
<link>
file          = links.txt
radius        = 0.98r
bezier_radius = 0.15r
thickness     = 2
ribbon        = yes
flat          = yes
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
    print(f"wrote {outdir}/circos.conf (+ ideogram.conf, ticks.conf)", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("hsp_detail_tsv")
    p.add_argument(
        "fasta_dir", help="Directory of per-genome FASTAs (e.g. vibriowatch_cholera_fasta/fasta_files_clean)"
    )
    p.add_argument(
        "--anchor",
        default="c6706",
        help="Anchor genome id, must appear as a 'genome' value in hsp_detail.tsv (default: c6706)",
    )
    p.add_argument(
        "--target-group",
        default="7PET",
        help="The group the anchor belongs to / markers are specific to (default: 7PET)",
    )
    p.add_argument(
        "--n-lineages",
        type=int,
        default=8,
        help="Max number of representative lineages to plot (default: 8, matches LINEAGE_COLORS length)",
    )
    p.add_argument(
        "--representatives",
        default=None,
        help="Pin exact genome IDs instead of auto-selecting, as group=genome,group=genome,...",
    )
    p.add_argument(
        "-o", "--outdir", default="circos_tier3_generated", help="Output directory (default: circos_tier3_generated)"
    )
    p.add_argument(
        "--etc-dir",
        default=None,
        help=(
            "Existing circos etc/ dir to reference from circos.conf "
            "(default: <outdir>/etc, copy it in yourself or pass the original circos_tier3/etc)"
        ),
    )
    p.add_argument("--image-name", default="circos_tier3", help="Output image basename (default: circos_tier3)")
    p.add_argument(
        "--render", action="store_true", help="Also run `circos -conf circos.conf` (needs the circos module loaded)"
    )
    args = p.parse_args()

    if args.n_lineages > len(LINEAGE_COLORS):
        sys.exit(
            f"--n-lineages {args.n_lineages} exceeds the fixed categorical palette length ({len(LINEAGE_COLORS)}); "
            "add more colours to LINEAGE_COLORS deliberately rather than cycling."
        )

    pinned = None
    if args.representatives:
        pinned = {}
        for pair in args.representatives.split(","):
            group, genome = pair.split("=", 1)
            pinned[group] = genome

    rows = load_hsp_detail(args.hsp_detail_tsv)
    if not rows:
        sys.exit(f"No rows in {args.hsp_detail_tsv}")

    representatives = pick_representatives(rows, args.anchor, args.target_group, args.n_lineages, pinned)
    if not representatives:
        sys.exit(
            "No representative lineages found -- check --anchor/--target-group match "
            "the table's 'genome'/'group' columns."
        )
    print(f"representatives: {representatives}", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genomes = [(args.anchor, args.anchor, ANCHOR_COLOR)]
    genomes += [(genome_id, group, OTHER_COLOR) for group, genome_id in representatives.items()]

    # Only keep contigs that actually carry a hit -- these are draft assemblies
    # that can run to hundreds of contigs; the original hand-built figure only
    # ever showed the handful with real data (see build_karyotype's docstring).
    relevant_genomes = {args.anchor} | set(representatives.values())
    used_contigs = defaultdict(set)
    for r in rows:
        if r["genome"] in relevant_genomes:
            used_contigs[r["genome"]].add(r["contig"])

    contig_to_chr = build_karyotype(genomes, args.fasta_dir, outdir, used_contigs)
    build_links(rows, args.anchor, representatives, contig_to_chr, outdir)

    etc_dir = args.etc_dir or (outdir / "etc")
    if not Path(etc_dir).exists():
        sys.exit(
            f"{etc_dir} not found -- pass --etc-dir pointing at an existing circos etc/ "
            "(e.g. the original circos_tier3/etc), or copy one in first."
        )
    build_conf(outdir, etc_dir, args.image_name)

    if args.render:
        subprocess.run(["circos", "-conf", "circos.conf"], cwd=outdir, check=True)
        print(f"rendered {outdir}/{args.image_name}.png / .svg", file=sys.stderr)
    else:
        print(f"run: cd {outdir} && circos -conf circos.conf", file=sys.stderr)


if __name__ == "__main__":
    main()
