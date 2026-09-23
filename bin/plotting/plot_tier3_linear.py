#!/usr/bin/env python3
"""Linear genome map of Tier 3 blastn hits.

Where each marker sits on the V. cholerae type-strain reference (N16961, 2 chromosomes)
and how far it leaks: one lollipop per marker at its genomic position, height / colour =
number of non-7PET genomes that also carry it (from the Tier 3 verdict table).

    python3 plot_tier3_linear.py <blast_hits.tsv> <tier3_validation.tsv> <out_dir>
    -> <out_dir>/panels/panel_E_genome_map.png / .pdf  (+ note in README.md)
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
                     "svg.fonttype": "none", "pdf.fonttype": 42})

REF = "n16961"                                     # type strain, El Tor O1 7PET
CHR = {"NZ_LT906614": ("chromosome 1", 2961182),   # sizes from the assembly
       "NZ_LT906615": ("chromosome 2", 1072331)}


def main(blast_path, val_path, out_dir):
    d = Path(out_dir) / "panels"
    d.mkdir(parents=True, exist_ok=True)

    # per-marker leak count from the Tier 3 verdict table
    leak = {}
    for r in csv.DictReader(open(val_path), delimiter="\t"):
        leak[r["marker_id"]] = int(r["n_nontarget_genomes"])

    # marker position on the reference (one HSP per marker)
    pos = {}   # marker -> (chrkey, midpoint)
    with open(blast_path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            sub = row[1]
            if sub.lower().startswith(REF + "__"):
                acc = sub.split("__", 1)[1].split(".")[0]
                if acc in CHR:
                    s, e = int(row[8]), int(row[9])
                    pos[row[0]] = (acc, (min(s, e) + max(s, e)) // 2)

    fig, axs = plt.subplots(2, 1, figsize=(14, 6.2), gridspec_kw=dict(hspace=0.55))
    vmax = max(leak.values()) if leak else 1
    norm = mcolors.SymLogNorm(linthresh=10, vmin=0, vmax=vmax)
    cmap = plt.get_cmap("YlOrRd")

    for ax, (acc, (title, size)) in zip(axs, CHR.items()):
        ms = [(p, leak.get(m, 0), m) for m, (c, p) in pos.items() if c == acc]
        ms.sort()
        for p, lk, m in ms:
            ax.vlines(p / 1e6, 0, max(lk, 0.8), color=cmap(norm(lk)), lw=2.2)
            ax.scatter([p / 1e6], [max(lk, 0.8)], s=42, color=cmap(norm(lk)),
                       edgecolor="#333", linewidth=0.5, zorder=3)
        ax.set_xlim(-0.03, size / 1e6 + 0.03)
        ax.set_ylim(0, vmax * 1.15 if vmax else 1)
        ax.set_yscale("symlog", linthresh=10)
        ax.set_title(f"{title}   ({size/1e6:.2f} Mb,  {len(ms)} markers)", fontsize=13, loc="left")
        ax.set_xlabel("position (Mb)", fontsize=12)
        ax.set_ylabel("non-7PET genomes\ncarrying the marker", fontsize=11)
        ax.tick_params(labelsize=11)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=axs, fraction=0.03, pad=0.02)
    cb.set_label("non-7PET genomes carrying the marker", fontsize=11)
    cb.ax.tick_params(labelsize=10)

    fig.savefig(d / "panel_E_genome_map.png", dpi=250, facecolor="white", bbox_inches="tight")
    fig.savefig(d / "panel_E_genome_map.pdf", facecolor="white", bbox_inches="tight")
    plt.close(fig)

    # append to README
    readme = d / "README.md"
    n_loci = len({(c, round(p, -3)) for m, (c, p) in pos.items()})
    extra = (f"\n- **panel_E_genome_map** — every marker's position on the N16961 type-strain "
             f"reference (2 chromosomes). One lollipop per marker; height & colour = number of "
             f"non-7PET genomes that also carry it. The {len(pos)} markers collapse to ~{n_loci} "
             f"loci — the vertical stacks are groups of near-duplicate candidate unitigs from one "
             f"underlying sequence. Tall red stacks are the most non-specific loci.\n")
    if readme.exists():
        readme.write_text(readme.read_text().rstrip() + "\n" + extra)
    else:
        readme.write_text("# Tier 3 panels\n" + extra)
    print(f"wrote panel_E_genome_map to {d}  ({len(pos)} markers on {REF})")


if __name__ == "__main__":
    main(*sys.argv[1:4])
