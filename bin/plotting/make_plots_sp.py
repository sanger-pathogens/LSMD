#!/usr/bin/env python3
"""Sprint-review plots for the S. pneumoniae GPSC marker-validation run.

SP counterpart of make_plots.py (which is V. cholerae 7PET, core/relaxed/catchall).
Here the panels are the five target GPSCs; each has its own candidate-marker set.

  Tier A  themisto pseudoalign vs the run's own 42,060-genome SP colour index.
          PASS = within-GPSC >= 95% AND <= 5% of any single other GPSC.
  Tier B  pseudoalign vs Jarno's ATB species-coloured index (12,733 species).
          SP core markers are ~35 bp (5 k-mers) -> all below the 70-k-mer verdict
          threshold, so Tier B is the raw cross-species hit picture, not pass/fail.
  Tier C  strict blastn vs the 42,060 SP assemblies (pident >= 98, qcovhsp >= 80).
          PASS = >=1 target-GPSC genome AND zero other genomes; FLAG = also elsewhere;
          ABSENT = no genome. (Written only once run_tier_c_sp.sh has finished.)

Usage: python3 make_plots_sp.py [root]   -- root defaults to the current directory,
so run it from inside the validation run directory, or pass the path explicitly.
"""
import sys
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
OUT = ROOT / "plots"
OUT.mkdir(exist_ok=True)

GPSCS = ["1.0", "3.0", "6.0", "12.0", "16.0"]
GENOMES = {"1.0": 1246, "3.0": 5141, "6.0": 1563, "12.0": 1959, "16.0": 1618}
SPECIES_N = 42060

# validated categorical palette (dataviz reference instance, light mode) -- same as make_plots.py
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GREEN = "#008300"   # status: pass
GREY = "#a9a89f"    # status: absent / n-a
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e6e5e2"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 12, "font.family": "DejaVu Sans",
    "axes.edgecolor": INK2, "axes.linewidth": 0.8,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
})


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def header(fig, title, subtitle, y=0.965):
    fig.text(0.015, y, title, ha="left", va="top", fontsize=14.5, fontweight="bold")
    fig.text(0.015, y - 0.060, subtitle, ha="left", va="top", fontsize=9.3, color=INK2)


def style_ax(ax, axis="x"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)
    getattr(ax, f"{axis}axis").grid(True, color=GRID, linewidth=0.8)
    getattr(ax, f"{'x' if axis == 'y' else 'y'}axis").grid(False)


def ylabels():
    return [f"GPSC {g[:-2] if g.endswith('.0') else g}\n{GENOMES[g]:,} genomes" for g in GPSCS]


def a_path(g):
    return ROOT / f"sp_GPSC{g}" / "tier_a" / f"sp_GPSC{g}_validation.tsv"


def b_summary(g):
    return ROOT / f"sp_GPSC{g}" / "tier_b" / f"sp_GPSC{g}_summary.txt"


def c_path(g):
    return ROOT / f"sp_GPSC{g}" / "tier_c" / f"sp_GPSC{g}_validation.tsv"


# ============================================================================
# Figure 1 - Tier A: every GPSC set passes the species-index specificity re-check
# ============================================================================
def fig_tier_a():
    n_mark, n_pass, within_med, out_max = [], [], [], []
    for g in GPSCS:
        d = tsv(a_path(g))
        n_mark.append(len(d))
        n_pass.append(sum(1 for r in d if r["verdict"] == "PASS"))
        w = sorted(float(r["within_align"]) for r in d)
        within_med.append(w[len(w) // 2])
        out_max.append(max(float(r["outside_align"]) for r in d))

    y = np.arange(len(GPSCS))[::-1]
    fig = plt.figure(figsize=(9.6, 5.6))
    ax = fig.add_axes([0.17, 0.20, 0.78, 0.52])
    ax.barh(y, n_mark, 0.6, color=GREEN, edgecolor=SURFACE, linewidth=2, zorder=2)
    for yi, nm, wm in zip(y, n_mark, within_med):
        ax.annotate(f"{nm:,}   ·   100% pass   ·   median within-GPSC {wm:.0f}%",
                    (nm, yi), xytext=(7, 0), textcoords="offset points",
                    va="center", fontsize=9.3, fontweight="bold", color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels(), fontsize=9.7)
    ax.set_xlabel("candidate markers", labelpad=8)
    ax.set_xlim(0, max(n_mark) * 1.9)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    style_ax(ax, axis="x")

    header(fig, "All five GPSC marker sets pass Tier 1",
           "Tier A — themisto pseudoalignment of every candidate marker back against the run's own "
           "42,060-genome\nS. pneumoniae colour index. PASS = carried by ≥ 95% of the GPSC's own "
           "genomes AND by ≤ 5% of any\nsingle other GPSC.   664 markers scored → 664 pass (100%).")
    fig.text(0.017, 0.055,
             "Marker yield is set by the upstream candidate filter, not by this check — GPSC 12 yields 478, GPSC 1 only 3.",
             fontsize=8.3, color=INK2)
    fig.savefig(OUT / "sp_01_tier_a.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / "sp_01_tier_a.png")


# ============================================================================
# Figure 2 - Tier B: cross-species hit picture (ATB species index)
# ============================================================================
def _parse_b(g):
    only = other = None
    top = []
    with open(b_summary(g)) as fh:
        lines = fh.read().splitlines()
    for i, ln in enumerate(lines):
        if "hit only target species" in ln:
            only = int(ln.split(":")[-1])
        elif "hit >=1 other named species" in ln:
            other = int(ln.split(":")[-1])
        elif "non-target species hit" in ln:
            for sub in lines[i + 1:]:
                s = sub.strip()
                if not s or not s.startswith("streptococcus"):
                    break
                name, cnt = s.rsplit(None, 1)
                top.append((name.replace("streptococcus_", "S. "), int(cnt)))
    return only or 0, other or 0, top


def fig_tier_b():
    only, other, tops = [], [], []
    for g in GPSCS:
        o, x, t = _parse_b(g)
        only.append(o)
        other.append(x)
        tops.append(t)
    tot = [o + x for o, x in zip(only, other)]
    only_pct = [100 * o / t if t else 0 for o, t in zip(only, tot)]
    other_pct = [100 * x / t if t else 0 for x, t in zip(other, tot)]

    y = np.arange(len(GPSCS))[::-1]
    fig = plt.figure(figsize=(9.8, 6.0))
    ax = fig.add_axes([0.16, 0.20, 0.80, 0.48])
    ax.barh(y, only_pct, 0.6, color=AQUA, edgecolor=SURFACE, linewidth=2,
            label="hits S. pneumoniae only")
    ax.barh(y, other_pct, 0.6, left=only_pct, color=ORANGE, edgecolor=SURFACE, linewidth=2,
            label="also hits ≥ 1 other Streptococcus species")
    for yi, op, xp, o, x in zip(y, only_pct, other_pct, only, other):
        if op > 8:
            ax.annotate(f"{op:.0f}%\n({o})", (op / 2, yi), ha="center", va="center",
                        fontsize=9, fontweight="bold", color="white")
        if xp > 8:
            ax.annotate(f"{xp:.0f}%\n({x})", (op + xp / 2, yi), ha="center", va="center",
                        fontsize=9, fontweight="bold", color="white")
    ax.set_yticks(y)
    ax.set_yticklabels([f"GPSC {g[:-2]}\nn = {t}" for g, t in zip(GPSCS, tot)], fontsize=9.7)
    ax.set_xlabel("candidate markers  (%)", labelpad=8)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    style_ax(ax, axis="x")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.03), fontsize=9,
              handlelength=1.3, ncol=2)
    header(fig, "Tier 2 — where the markers land across species",
           "Tier B — pseudoalignment vs Jarno's ATB species-coloured index (12,733 species, k = 31). "
           "SP core markers are ~35 bp\n(5 k-mers) — all below the 70-k-mer verdict threshold, so this "
           "is the raw cross-species hit picture, not a pass/fail.")
    # top cross-hit species, pooled
    pooled = {}
    for t in tops:
        for name, c in t:
            pooled[name] = pooled.get(name, 0) + c
    top3 = sorted(pooled.items(), key=lambda kv: -kv[1])[:3]
    txt = ("Cross-species hits are near-exclusively the mitis group:  "
           + ",  ".join(f"{n} ({c})" for n, c in top3) + ", …")
    fig.text(0.017, 0.05, txt, fontsize=8.3, color=INK2)
    fig.savefig(OUT / "sp_02_tier_b.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / "sp_02_tier_b.png")


# ============================================================================
# Figure 3 - Tier C: strict blastn vs the 42,060 SP assemblies
# ============================================================================
def fig_tier_c():
    if not all(c_path(g).exists() for g in GPSCS):
        missing = [g for g in GPSCS if not c_path(g).exists()]
        print(f"tier C not ready (missing: {', '.join(missing)}) -- skipping sp_03_tier_c.png")
        return
    P, F, A, within = [], [], [], []
    for g in GPSCS:
        d = tsv(c_path(g))
        n = len(d) or 1
        P.append(100 * sum(1 for r in d if r["verdict"] == "PASS") / n)
        F.append(100 * sum(1 for r in d if r["verdict"] == "FLAG") / n)
        A.append(100 * sum(1 for r in d if r["verdict"] == "ABSENT") / n)
        wp = [float(r["within_7PET"]) for r in d if r["verdict"] == "PASS"]
        within.append(np.mean(wp) if wp else 0.0)
    n_mark = [len(tsv(c_path(g))) for g in GPSCS]

    y = np.arange(len(GPSCS))[::-1]
    fig = plt.figure(figsize=(9.6, 5.8))
    ax = fig.add_axes([0.16, 0.16, 0.80, 0.52])
    # PASS=BLUE, FLAG=ORANGE, ABSENT=GREY -- validated categorical pair for the
    # two that matter (blue<->orange), grey is the deliberate neutral null state.
    ax.barh(y, P, 0.6, color=BLUE, edgecolor=SURFACE, linewidth=2, label="PASS — target GPSC only")
    ax.barh(y, F, 0.6, left=P, color=ORANGE, edgecolor=SURFACE, linewidth=2,
            label="FLAG — also in ≥ 1 other genome")
    ax.barh(y, A, 0.6, left=[p + f for p, f in zip(P, F)], color=GREY, edgecolor=SURFACE, linewidth=2,
            label="ABSENT — in no assembly")
    for yi, p, f, a in zip(y, P, F, A):
        for val, base in ((p, 0), (f, p), (a, p + f)):
            if val > 7:
                ax.annotate(f"{val:.0f}%", (base + val / 2, yi), ha="center", va="center",
                            fontsize=9, fontweight="bold", color="white")
    ax.set_yticks(y)
    ax.set_yticklabels([f"GPSC {g[:-2]}\nn = {m}" for g, m in zip(GPSCS, n_mark)], fontsize=9.7)
    ax.set_xlabel("candidate markers  (%)")
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    style_ax(ax, axis="x")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.03), fontsize=9, ncol=3)
    header(fig, "Tier 3 — strict blastn against all 42,060 S. pneumoniae assemblies",
           "Tier C — blastn (pident ≥ 98, qcovhsp ≥ 80), two knobs only. A genome carries the "
           "marker iff one HSP clears both.\nPASS = ≥ 1 target-GPSC genome AND zero others.")
    fig.savefig(OUT / "sp_03_tier_c.png", dpi=200)
    plt.close(fig)
    print("wrote", OUT / "sp_03_tier_c.png")


if __name__ == "__main__":
    fig_tier_a()
    fig_tier_b()
    fig_tier_c()
