"""Figures for the embedding supplementary (notebook 08).

One module per output family, under the package rather than beside the notebook,
so there is a single import root. Each function RETURNS its figure; saving
belongs to `biokgsuite.io.Outputs`, as in `viz.dimension04` and `viz.dimension05`.

All of these are driven by cached CSV tables written during the analysis, so
they regenerate without retraining any embeddings.

Colour follows the entity. A graph is always its ``KG_PALETTE`` colour. The
trained scorers (TransE, RotatE) are models, not graphs, so they take the
model palette shared with the grounding figures, each with its own marker
shape; the text-only name prior (Gemma) is the baseline the other two are
measured against and is drawn in ``NEUTRAL`` grey like every other baseline.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from biokgsuite.viz.grounding import MODEL_MARKERS, MODEL_PALETTE
from biokgsuite.viz.primitives import facet_grid, panel_letters, ref_line
from biokgsuite.viz.style import (DOUBLE_COL_W, FAINT, KG_NAMES, KG_ORDER, KG_PALETTE, MUTED,
                                  NEUTRAL, SINGLE_COL_W, TEXT_COLOR, clean_ax, facet_title,
                                  fig_size, setup_style)

setup_style()

MCOLORS = {'TransE': MODEL_PALETTE[0], 'RotatE': MODEL_PALETTE[1], 'Gemma': NEUTRAL}
MMARKERS = {'TransE': MODEL_MARKERS[0], 'RotatE': MODEL_MARKERS[1], 'Gemma': MODEL_MARKERS[2]}
KG_LABEL = dict(KG_NAMES)
KG_COLORS = {k: KG_PALETTE.get(k, NEUTRAL) for k in KG_LABEL}
NQ = {'primekg': 'readable', 'hetionet': 'readable', 'drkg': 'id-only',
      'openbilink': 'id-only', 'biokg': 'mixed names', 'matrix': 'mixed'}

BAR_W = 0.6
ERR_KW = dict(elinewidth=0.75, capsize=0, zorder=3)


def _kg_order(kgs):
    """Graphs in the fixed order of the paper; unknown keys follow, as given."""
    kgs = list(kgs)
    return [k for k in KG_ORDER if k in kgs] + [k for k in kgs if k not in KG_ORDER]


def _grouped_bars(ax, kgs, series, colours):
    """Bars per graph, one per series, the group 0.6 of the slot."""
    x = np.arange(len(kgs))
    w = BAR_W / max(len(series), 1)
    for i, (name, means, errs) in enumerate(series):
        ax.bar(x + (i - (len(series) - 1) / 2) * w, means, w, color=colours[name], label=name,
               zorder=2, yerr=errs, error_kw=dict(ecolor=TEXT_COLOR, **ERR_KW))
    ax.set_xticks(x)
    ax.set_xlim(-0.6, len(kgs) - 0.4)
    ax.set_xticklabels([KG_LABEL.get(k, k) for k in kgs])
    ax.tick_params(axis='x', length=0)


def _name_quality(ax, kgs):
    """The naming regime of each graph, in 6 pt muted italics under its name:
    a text-only prior can only score what it can read."""
    for i, k in enumerate(kgs):
        if k in NQ:
            ax.annotate(NQ[k], (i, 0), xycoords=('data', 'axes fraction'), xytext=(0, -14),
                        textcoords='offset points', ha='center', va='top', fontsize=6,
                        color=MUTED, style='italic', annotation_clip=False)


# ───────────────────────── 08 · resampled AUROC (grouped bars) ──────────────────
def fig08_resampled_auroc(resampled_csv):
    STRAT = 'type-constrained'
    MODELS = ['TransE', 'RotatE', 'Gemma']
    df = pd.read_csv(resampled_csv)
    df = df[df.strategy == STRAT]
    models = [m for m in MODELS if m in set(df.model)]

    def _stats(kg, m):
        v = df[(df.kg == kg) & (df.model == m)]['auroc'].dropna()
        se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
        return v.mean(), 1.96 * se

    kgs = _kg_order(df.kg.unique())
    fig, ax = plt.subplots(figsize=fig_size(1 if len(kgs) <= 3 else 2, 1), constrained_layout=True)
    series = [(m, [_stats(k, m)[0] for k in kgs], [_stats(k, m)[1] for k in kgs]) for m in models]
    _grouped_bars(ax, kgs, series, MCOLORS)
    _name_quality(ax, kgs)
    ref_line(ax, 0.5, 'chance')
    ax.set_ylim(0, 1.0)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    clean_ax(ax, ylabel='AUROC, type-constrained negatives')
    ax.legend(loc='upper left', ncol=len(models))
    return fig


# ───────────────────────── 08 · AUROC lift over Gemma prior ─────────────────────
def fig08_lift_over_gemma(resampled_csv):
    TRAINED = ['TransE', 'RotatE']
    BASELINE = 'Gemma'
    df = pd.read_csv(resampled_csv)
    df = df[df.strategy == 'type-constrained']
    trained = [m for m in TRAINED if m in set(df.model)]
    rep_col = 'rerun' if 'rerun' in df.columns else None

    def _lift(kg, model):
        v = df[(df.kg == kg) & (df.model == model)]['auroc'].dropna()
        g = df[(df.kg == kg) & (df.model == BASELINE)]['auroc'].dropna()
        if len(v) == 0 or len(g) == 0:
            return np.nan, 0.0
        if rep_col:
            a = df[(df.kg == kg) & (df.model == model)][[rep_col, 'auroc']]
            b = df[(df.kg == kg) & (df.model == BASELINE)][[rep_col, 'auroc']]
            dd = a.merge(b, on=rep_col, suffixes=('_m', '_g'))
            dd = (dd['auroc_m'] - dd['auroc_g']).dropna()
            if len(dd) > 1:
                return dd.mean(), 1.96 * dd.std(ddof=1) / np.sqrt(len(dd))
        se_v = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
        se_g = g.std(ddof=1) / np.sqrt(len(g)) if len(g) > 1 else 0.0
        return v.mean() - g.mean(), 1.96 * np.sqrt(se_v ** 2 + se_g ** 2)

    kgs = _kg_order([k for k in df.kg.unique()
                     if not df[(df.kg == k) & (df.model == BASELINE)].empty
                     and any(not df[(df.kg == k) & (df.model == m)].empty for m in trained)])
    fig, ax = plt.subplots(figsize=fig_size(1 if len(kgs) <= 3 else 2, 1), constrained_layout=True)
    series = [(m, [_lift(k, m)[0] for k in kgs], [_lift(k, m)[1] for k in kgs]) for m in trained]
    _grouped_bars(ax, kgs, series, MCOLORS)
    lows = [mu - e for _, means, errs in series for mu, e in zip(means, errs)]
    highs = [mu + e for _, means, errs in series for mu, e in zip(means, errs)]
    ymin, ymax = min([0] + lows), max([0] + highs)
    pad = 0.08 * (ymax - ymin or 1)
    # Head room above the tallest bar for the legend.
    ax.set_ylim(ymin - pad, ymax + 3 * pad)
    ref_line(ax, 0, 'Gemma prior')
    clean_ax(ax, ylabel='AUROC gain over Gemma prior')
    ax.legend(loc='upper left', ncol=len(trained))
    return fig


# ──────────────── 08 · heuristic vs embedding AUROC scatter (per strategy) ───────
def fig08_heuristic_scatter(single_csv, resampled_csv):
    sdf = pd.read_csv(single_csv)
    heur = sdf.groupby(["kg", "strategy"])["heuristic_auroc"].first().reset_index()
    rdf = pd.read_csv(resampled_csv)
    emb = (rdf[rdf.model.isin(["TransE", "RotatE"])]
           .groupby(["kg", "model", "strategy"])["auroc"]
           .agg(emb_auroc="mean",
                emb_lo=lambda x: np.percentile(x, 2.5),
                emb_hi=lambda x: np.percentile(x, 97.5),
                n="count")
           .reset_index())
    comp_df = (emb.merge(heur, on=["kg", "strategy"], how="inner")
               .dropna(subset=["emb_auroc", "heuristic_auroc"]))

    STRATEGIES = [s for s in ["random", "degree-matched", "type-constrained"]
                  if s in set(comp_df.strategy)] or sorted(comp_df.strategy.unique())

    rows = []
    for m in ["TransE", "RotatE"]:
        for s in STRATEGIES:
            d = comp_df[(comp_df.model == m) & (comp_df.strategy == s)]
            if len(d) >= 3:
                rho, p = spearmanr(d.heuristic_auroc, d.emb_auroc)
                rows.append((m, s, rho, p))
    spearman_df = pd.DataFrame(rows, columns=["model", "strategy", "spearman_rho", "p_value"])

    vals = np.r_[comp_df.heuristic_auroc.values, comp_df.emb_auroc.values,
                 comp_df.emb_lo.values, comp_df.emb_hi.values]
    pad = 0.03 * (vals.max() - vals.min() or 1)
    LO, HI = vals.min() - pad, vals.max() + pad

    n = len(STRATEGIES)
    # Square panels: three across the double column are 55 mm on a side.
    row_h = (SINGLE_COL_W if n == 1 else DOUBLE_COL_W / min(n, 3)) + 0.25
    fig, axes, ncol = facet_grid(n, row_h=row_h)
    for ax_i, strat in enumerate(STRATEGIES):
        ax = axes[ax_i]
        ax.plot([LO, HI], [LO, HI], ls='--', color=FAINT, lw=0.6, zorder=1)
        for model_name in ("TransE", "RotatE"):
            sub = comp_df[(comp_df.model == model_name) & (comp_df.strategy == strat)].dropna(
                subset=["emb_auroc", "heuristic_auroc"])
            for _, r in sub.iterrows():
                ax.errorbar(r.heuristic_auroc, r.emb_auroc,
                            yerr=[[r.emb_auroc - r.emb_lo], [r.emb_hi - r.emb_auroc]],
                            fmt="none", ecolor=KG_COLORS.get(r.kg, NEUTRAL), elinewidth=0.75,
                            capsize=0, zorder=4)
            ax.scatter(sub.heuristic_auroc, sub.emb_auroc,
                       c=[KG_COLORS.get(k, NEUTRAL) for k in sub.kg],
                       marker=MMARKERS[model_name], s=14, zorder=5, edgecolors="white",
                       linewidths=0.4)
        _sr = spearman_df[(spearman_df.model == "TransE") & (spearman_df.strategy == strat)]
        if len(_sr):
            r0 = _sr.iloc[0]
            ax.text(0.04, 0.96, f"TransE ρ = {r0.spearman_rho:.2f}, P = {r0.p_value:.3f}",
                    transform=ax.transAxes, fontsize=6, va="top", color=MUTED)
        ax.set_xlim(LO, HI)
        ax.set_ylim(LO, HI)
        ax.set_aspect("equal", adjustable="box")
        clean_ax(ax, xlabel="Heuristic AUROC" if ax_i >= n - ncol else '',
                 ylabel="Embedding AUROC" if ax_i % ncol == 0 else '')
        if ax_i % ncol:
            ax.tick_params(labelleft=False)
        facet_title(ax, strat[:1].upper() + strat[1:])
    # Shape says which scorer, colour says which graph; both legends are
    # frameless and the graph legend lists only the graphs that were scored.
    axes[n - 1].legend(
        handles=[Line2D([0], [0], marker=MMARKERS[m], color=TEXT_COLOR, lw=0, ms=3.5,
                        label=m, markeredgecolor="white", markeredgewidth=0.4)
                 for m in ("TransE", "RotatE")],
        loc="lower right")
    present = _kg_order(comp_df.kg.unique())
    kg_handles = [Line2D([0], [0], marker="o", color=KG_COLORS.get(k, NEUTRAL), lw=0, ms=3.5,
                         markeredgecolor="white", markeredgewidth=0.4, label=KG_LABEL.get(k, k))
                  for k in present]
    fig.legend(handles=kg_handles, ncol=len(kg_handles), loc="outside lower center")
    if n > 1:
        panel_letters(fig, axes)
    return fig
