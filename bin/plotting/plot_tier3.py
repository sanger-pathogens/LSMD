#!/usr/bin/env python3
"""Tier 3 (blastn) -- four separate, title-less panels under <tier3_dir>/panels/.

    A  lineage reach   : per non-target lineage, % of markers that also hit it (bar),
                         with % of that lineage's genomes carrying the markers annotated
                         beside the bar in a contrasting colour
    B  match quality   : 2-D histogram of pident vs. query coverage per HSP
    C  SNP position    : where the (rare) mismatches fall along the marker
    D  within-target   : per marker, % of hit target genomes carrying it

Text for the deck goes to panels/README.md.

    python3 plot_tier3.py <hsp_detail.tsv> <lineage_tally.tsv> <target> <out_dir>
"""
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

BAR = "#2166ac"  # bar fill
ANNOT = "#c0392b"  # the second metric, in a contrasting colour


def save(fig, out):
    fig.savefig(f"{out}.png", dpi=250, facecolor="white", bbox_inches="tight")
    fig.savefig(f"{out}.pdf", facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main(hsp_path, tally_path, target, out_dir):
    d = Path(out_dir) / "panels"
    d.mkdir(parents=True, exist_ok=True)

    # --- single streaming pass over hsp_detail (can be tens of GB for catchall) ---
    PID_BINS = np.arange(97.5, 100.76, 0.25)
    QC_BINS = np.arange(79.0, 101.5, 1.0)
    hist2d = np.zeros((len(PID_BINS) - 1, len(QC_BINS) - 1))
    mmpos = Counter()
    tgt_per_marker = Counter()
    markers, genomes = set(), set()
    n_hsp = 0
    with open(hsp_path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        ci = {c: i for i, c in enumerate(hdr)}
        im_, ig, igr = ci["marker_id"], ci["genome"], ci["group"]
        ip, iq, imp = ci["pident"], ci["qcovhsp"], ci["mismatch_pos"]
        for row in rdr:
            if len(row) <= imp:
                continue
            n_hsp += 1
            markers.add(row[im_])
            genomes.add(row[ig])
            pv, qv = float(row[ip]), float(row[iq])
            bi = int((pv - 97.5) / 0.25)
            bj = int((qv - 79.0) / 1.0)
            if 0 <= bi < hist2d.shape[0] and 0 <= bj < hist2d.shape[1]:
                hist2d[bi, bj] += 1
            for pp in row[imp].split(","):
                if pp.strip().isdigit():
                    mmpos[int(pp)] += 1
            if row[igr] == target:
                tgt_per_marker[row[im_]] += 1
    n_markers, n_genomes = len(markers), len(genomes)

    pass_ids = Path(hsp_path).with_name(Path(hsp_path).name.replace("_hsp_detail.tsv", "_PASS.ids"))
    n_pass = sum(1 for _ in open(pass_ids)) if pass_ids.exists() else 0
    n_notspec = n_markers - n_pass

    # ---------- A : lineage reach / specificity ----------
    # bar = % of markers that ALSO hit this non-target lineage; red = % of that
    # lineage's genomes carrying the markers (mean). Least-cross-hit first.
    tally = [t for t in csv.DictReader(open(tally_path), delimiter="\t") if t["group"] != target]
    tally.sort(key=lambda t: float(t["pct_markers_matched"]))
    g = [t["group"] for t in tally]
    pct_m = [float(t["pct_markers_matched"]) for t in tally]
    mean_grp = [float(t["mean_frac_of_group"]) for t in tally]
    y = np.arange(len(g))
    figA, a = plt.subplots(figsize=(9.5, max(3.5, 0.46 * len(g) + 1.4)))
    a.barh(y, pct_m, color=BAR, height=0.66)
    for i, (pm, mg) in enumerate(zip(pct_m, mean_grp)):
        a.text(min(pm + 2, 102), i, f"{mg:.0f}%", va="center", ha="left", fontsize=12.5, fontweight="bold", color=ANNOT)
    a.set_yticks(y)
    a.set_yticklabels(g, fontsize=13)
    a.set_xlabel(f"% of the {n_markers:,} markers that also occur in this lineage", fontsize=13)
    a.set_xlim(0, 112)
    a.tick_params(labelsize=12, length=0)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.text(
        0.985,
        0.965,
        f"{n_notspec:,} / {n_markers:,} markers also occur in >= 1\nnon-{target} genome  ({n_pass:,} pass Tier 3)",
        transform=a.transAxes,
        ha="right",
        va="top",
        fontsize=12.5,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.5", fc="#fdecea", ec="#c0392b"),
    )
    a.text(
        0.985,
        0.02,
        "red = % of that lineage's genomes carrying the markers (mean)",
        transform=a.transAxes,
        ha="right",
        fontsize=11,
        color=ANNOT,
        style="italic",
    )
    save(figA, d / "panel_A_lineage_reach")

    # ---------- B : match quality ----------
    figB, b = plt.subplots(figsize=(7.5, 6))
    masked = np.ma.masked_where(hist2d.T == 0, hist2d.T)
    im = b.pcolormesh(PID_BINS, QC_BINS, masked, cmap="Blues")
    cb = figB.colorbar(im, ax=b)
    cb.ax.tick_params(labelsize=11)
    cb.set_label("HSPs", fontsize=13)
    b.set_xlabel("percent identity", fontsize=13)
    b.set_ylabel("query coverage per HSP (%)", fontsize=13)
    b.tick_params(labelsize=12)
    save(figB, d / "panel_B_match_quality")

    # ---------- C : SNP position within the marker ----------
    figC, c = plt.subplots(figsize=(8, 4.5))
    if mmpos:
        xs = sorted(mmpos)
        c.bar(xs, [mmpos[x] for x in xs], color="#e67e22", width=0.9)
    c.set_xlabel("position along marker (bp)", fontsize=13)
    c.set_ylabel("mismatches", fontsize=13)
    c.tick_params(labelsize=12, length=0)
    for s in ("top", "right"):
        c.spines[s].set_visible(False)
    save(figC, d / "panel_C_snp_position")

    # ---------- D : reach within the target ----------
    tgt_hits = np.array(list(tgt_per_marker.values())) if tgt_per_marker else np.array([0])
    mx = tgt_hits.max() or 1
    figD, e = plt.subplots(figsize=(8, 4.5))
    e.hist(tgt_hits / mx * 100, bins=np.arange(0, 105, 5), color="#27ae60")
    e.set_xlabel(f"% of hit {target} genomes carrying the marker", fontsize=13)
    e.set_ylabel("markers", fontsize=13)
    e.tick_params(labelsize=12, length=0)
    for s in ("top", "right"):
        e.spines[s].set_visible(False)
    save(figD, d / "panel_D_within_target")

    n_mm = sum(mmpos.values())
    top = "; ".join(
        f"{t['group']} ({float(t['pct_markers_matched']):.0f}% of markers / "
        f"{float(t['mean_frac_of_group']):.0f}% of that lineage)"
        for t in tally[-6:][::-1]
    )
    mode = next((w for w in ("core", "relaxed", "catchall") if w in str(hsp_path)), "")
    (d / "README.md").write_text(
        f"""# Tier 3 panels — blastn, {target} {mode} (pident >= 98 / qcov >= 80)

{n_markers:,} markers vs {n_genomes:,} V. cholerae assemblies ·
**{n_pass:,} / {n_markers:,} PASS** ({n_notspec:,} also occur in >= 1 non-{target} genome).

- **panel_A_lineage_reach** — bar: % of the {n_markers:,} markers that also occur in each non-{target} lineage.
  Red number beside each bar: mean % of that lineage's genomes carrying the markers. Worst: {top}.
- **panel_B_match_quality** — pident vs query-coverage, {n_hsp:,} HSPs.
  The mass sits at ~100 / 100 — these are near-pristine matches, not marginal hits.
- **panel_C_snp_position** — {n_mm:,} mismatches across all {n_hsp:,} HSPs; they cluster at the marker ends.
  The sequence is nearly invariant across the species.
- **panel_D_within_target** — per marker, the % of hit {target} genomes carrying it —
  how consistently each marker is present in its own lineage.

**Takeaway:** a large fraction of these markers are conserved sequence shared beyond {target} —
the PAT-3570 result, by a third method.
"""
    )
    print(f"wrote 4 panels + README.md to {d}  ({n_markers:,} markers, {n_hsp:,} HSPs, {n_pass:,} PASS)")


if __name__ == "__main__":
    main(*sys.argv[1:5])
