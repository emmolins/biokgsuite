#!/usr/bin/env python3
"""
Benchmark construction figure for the cleaned (v5) evaluation set.

Replaces the earlier 253 -> 234 -> 143 -> 116 figure. Two changes of substance:
the -91/-27 split is gone because it could not be reproduced from the stored
tables (see METHODS_BENCHMARK_CONSTRUCTION.md section 9), and the funnel now
continues past 116 through deduplication and register verification to 80.
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK, MUTED, FAINT = "#1a1a1a", "#5b5b5b", "#8a8a8a"
DROP = "#c0392b"
KG = {"PrimeKG": "#1f6feb", "DRKG": "#3fa34d", "BioKG": "#8250df"}
SURF, RULE = "#fcfcfb", "#dcdcd8"

STAGES = [
    ("253", "candidate pairs",
     "new-indication approvals 2023–2026, four agencies", False),
    ("234", "harness-ready", "resolved drug and disease identifiers", False),
    ("116", "evaluable", "as previously reported", False),
    ("99", "unique drug–disease pairs", "one query per pair", False),
    ("80", "verified evaluation set",
     "every row cited to a regulator record", True),
]
# Wrapped by hand: text drawn in panel a data coordinates is not clipped, so an
# over-long line silently runs across the right-hand panels.
DROPS = [
    ("−19", "unresolved IDs, or\nplaceholder indication text"),
    ("−118", "answer leakage, or drug/disease\nabsent from all three KGs"),
    ("−17", "duplicate pair: the same\napproval recorded twice"),
    ("−19", "5 contradicted by the register\n4 not an indication change\n10 population or severity expansion"),
]


def panel_a(ax):
    ax.set_xlim(0, 11.4); ax.set_ylim(0, 10.6); ax.axis("off")
    top, h, gap = 9.35, 1.28, 0.72
    for i, (n, title, sub, hero) in enumerate(STAGES):
        y = top - i * (h + gap)
        ax.add_patch(FancyBboxPatch(
            (0, y - h), 6.25, h, boxstyle="round,pad=0.02,rounding_size=0.06",
            linewidth=2.0 if hero else 1.0,
            edgecolor=INK if hero else RULE,
            facecolor="#f4f4f1" if hero else SURF, zorder=2))
        ax.text(0.34, y - 0.52, n, fontsize=19 if hero else 17,
                fontweight="bold", color=INK, va="center", zorder=3)
        ax.text(1.62, y - 0.42, title, fontsize=9.5,
                fontweight="bold" if hero else "normal", color=INK,
                va="center", zorder=3)
        ax.text(1.62, y - 0.86, sub, fontsize=7.8, color=MUTED, va="center", zorder=3)
        if i < len(STAGES) - 1:
            ymid = y - h - gap / 2
            ax.add_patch(FancyArrowPatch(
                (0.9, y - h - 0.05), (0.9, y - h - gap + 0.05),
                arrowstyle="-|>", mutation_scale=9, linewidth=1.0,
                color=FAINT, zorder=2))
            dn, dl = DROPS[i]
            ax.plot([1.15, 6.15], [ymid, ymid], linestyle=(0, (2, 2)),
                    linewidth=0.9, color=DROP, zorder=1)
            ax.text(6.4, ymid, dn, fontsize=9.5, fontweight="bold",
                    color=DROP, va="center", zorder=3)
            ax.text(7.35, ymid, dl, fontsize=7.3, color=MUTED,
                    va="center", linespacing=1.5, zorder=3)


def panel_b(ax, cov, allthree, n_map):
    names = list(cov.keys()); vals = [cov[k] for k in names]
    y = range(len(names))[::-1]
    ax.barh(list(y), [n_map] * len(names), height=0.44,
            color="#ececea", edgecolor="none", zorder=1)
    for yy, nm, v in zip(y, names, vals):
        ax.barh(yy, v, height=0.44, color=KG[nm], edgecolor="none", zorder=2)
        ax.text(-1.4, yy, nm, ha="right", va="center", fontsize=9, color=INK)
        ax.text(v + 1.4, yy, str(v), va="center", fontsize=9,
                fontweight="bold", color=INK)
    ax.set_xlim(0, n_map + 6); ax.set_ylim(-1.5, len(names) - 0.3)
    ax.axis("off")
    ax.text(0, -0.70, f"of {n_map} pairs resolvable in the coverage audit; "
                      f"7 combination products excluded",
            fontsize=7.3, color=MUTED)
    ax.text(0, -1.22, f"{allthree}", fontsize=15, fontweight="bold",
            color=INK, va="center")
    ax.text(5.2, -1.22, "covered by all three: the common support on which\n"
                        "evidence quality is comparable across graphs",
            fontsize=7.3, color=MUTED, linespacing=1.5, va="center")


def stacked(ax, yy, parts, colors, total, label):
    x = 0
    for (nm, v), c in zip(parts, colors):
        ax.barh(yy, v, left=x, height=0.4, color=c, edgecolor=SURF,
                linewidth=1.6, zorder=2)
        if v / total > 0.1:
            ax.text(x + v / 2, yy, str(v), ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white", zorder=3)
        x += v
    ax.text(-1.2, yy, label, ha="right", va="center", fontsize=8.5, color=INK)


def panel_c(ax, ev, ver):
    stacked(ax, 1, ev, ["#1f6feb", "#3fa34d", "#8250df"], 80, "register")
    stacked(ax, 0, ver, ["#1a1a1a", "#b9b9b4"], 80, "verification")
    ax.set_xlim(0, 86); ax.set_ylim(-0.95, 1.5); ax.axis("off")
    leg = [("Drugs@FDA", "#1f6feb"), ("EMA EPAR", "#3fa34d"),
           ("PMDA list", "#8250df"), ("disease-level", "#1a1a1a"),
           ("drug + window", "#b9b9b4")]
    x = 0
    for nm, c in leg:
        ax.add_patch(plt.Rectangle((x, -0.62), 1.7, 0.15, color=c, zorder=3))
        ax.text(x + 2.3, -0.54, nm, fontsize=7.1, color=MUTED, va="center")
        x += 2.3 + len(nm) * 1.02 + 2.6


def main():
    cj = json.loads(Path("results/tables/v5_strict_coverage.json").read_text())
    cov, allthree, n_map = cj["coverage"], cj["all_three"], cj["n_mapped"]

    fig = plt.figure(figsize=(13.2, 6.0), dpi=200)
    fig.patch.set_facecolor(SURF)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.30, 0.92],
                          height_ratios=[1.0, 0.52],
                          left=0.030, right=0.985, top=0.905, bottom=0.085,
                          wspace=0.20, hspace=1.05)
    axa = fig.add_subplot(gs[:, 0]); axa.set_facecolor(SURF)
    axb = fig.add_subplot(gs[0, 1]); axb.set_facecolor(SURF)
    axc = fig.add_subplot(gs[1, 1]); axc.set_facecolor(SURF)

    # Panel letters and titles in FIGURE coordinates: placing them per-axes in
    # transAxes drifts with each panel's own label gutter and was what pushed
    # panel c's title into panel a's annotation column.
    for x, letter, title in ((0.030, "a", "Benchmark construction"),
                             (0.605, "b", "Knowledge-graph coverage of the 80"),
                             (0.605, "c", "Evidence behind each pair")):
        yy = 0.945 if letter == "a" else (0.945 if letter == "b" else 0.335)
        fig.text(x, yy, letter, fontsize=13, fontweight="bold", color=INK)
        fig.text(x + 0.022, yy, title, fontsize=10.5, fontweight="bold", color=INK)

    panel_a(axa)
    panel_b(axb, cov, allthree, n_map)
    panel_c(axc, [("Drugs@FDA", 47), ("EMA EPAR", 23), ("PMDA list", 10)],
            [("disease-level", 29), ("drug + window", 51)])

    out = Path("results/figures")
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"benchmark_construction_v5.{ext}",
                    facecolor=SURF, bbox_inches="tight")
    print("wrote results/figures/benchmark_construction_v5.png / .pdf")


if __name__ == "__main__":
    main()
