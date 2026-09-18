"""Figures for the embedding supplementary (notebook 08).

One module per output family, under the package rather than beside the notebook,
so there is a single import root. Each function RETURNS its figure; saving
belongs to `biokgsuite.io.Outputs`, as in `viz.dimension04` and `viz.dimension05`.

All of these are driven by cached CSV tables written during the analysis, so
they regenerate without retraining any embeddings.

Model-identity colours (TransE / RotatE / Gemma) are intentionally NOT the KG
palette: they identify methods, not graphs. KG colours in the scatter use the
canonical KG_PALETTE.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

from biokgsuite.viz.style import KG_PALETTE


MCOLORS = {'TransE': '#1f4e79', 'RotatE': '#6699cc', 'Gemma': '#999999'}
KG_LABEL = {'primekg': 'PrimeKG', 'hetionet': 'Hetionet', 'drkg': 'DRKG',
            'openbilink': 'OpenBioLink', 'biokg': 'BioKG', 'matrix': 'Matrix'}
KG_COLORS = {k: KG_PALETTE.get(k, '#888888') for k in KG_LABEL}
NQ = {'primekg': 'readable', 'hetionet': 'readable', 'drkg': 'id-only',
      'openbilink': 'id-only', 'biokg': 'mixed names', 'matrix': 'mixed'}


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

    kgs = sorted(df.kg.unique(), key=lambda k: max(_stats(k, m)[0] for m in models), reverse=True)
    x, w = np.arange(len(kgs)), 0.78 / len(models)

    fig, ax = plt.subplots(figsize=(7.0, 3.3))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color='#E8E8E8', lw=0.6, zorder=0)
    ax.xaxis.grid(False)
    for i, m in enumerate(models):
        means = [_stats(k, m)[0] for k in kgs]
        errs = [_stats(k, m)[1] for k in kgs]
        ax.bar(x + (i - (len(models) - 1) / 2) * w, means, w,
               color=MCOLORS[m], label=m, edgecolor='white', linewidth=0.4, zorder=2,
               yerr=errs, error_kw=dict(ecolor='#4a4a4a', elinewidth=0.8, capsize=2, capthick=0.8, zorder=3))
    ax.axhline(0.5, color='#BDBDBD', lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.annotate('chance', xy=(1.0, 0.5), xycoords=('axes fraction', 'data'),
                xytext=(4, 0), textcoords='offset points', ha='left', va='center',
                fontsize=6.5, color='#9a9a9a', annotation_clip=False)
    ax.set_ylim(0, 1.0)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_ylabel('AUROC (type-constrained negatives)')
    ax.set_xticks(x)
    ax.set_xlim(-0.6, len(kgs) - 0.4)
    ax.tick_params(bottom=False)
    ax.set_xticklabels([KG_LABEL[k] for k in kgs])
    for i, k in enumerate(kgs):
        ax.text(i, -0.10, NQ[k], transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=6.5, color='#8C8C8C', style='italic')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, ncol=len(models), loc='lower center', bbox_to_anchor=(0.5, 1.01),
              handlelength=1.1, handletextpad=0.5, columnspacing=1.8, borderpad=0)
    ax.set_title('Trained KG embeddings vs. a text-only name prior (Gemma)', loc='left', pad=20)
    fig.tight_layout()
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

    kgs = [k for k in df.kg.unique()
           if not df[(df.kg == k) & (df.model == BASELINE)].empty
           and any(not df[(df.kg == k) & (df.model == m)].empty for m in trained)]
    kgs.sort(key=lambda k: np.nanmean([_lift(k, m)[0] for m in trained]), reverse=True)
    x, w = np.arange(len(kgs)), 0.78 / len(trained)

    fig, ax = plt.subplots(figsize=(7.0, 3.3))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color='#E8E8E8', lw=0.6, zorder=0)
    ax.xaxis.grid(False)
    lows, highs = [], []
    for i, m in enumerate(trained):
        means = [_lift(k, m)[0] for k in kgs]
        errs = [_lift(k, m)[1] for k in kgs]
        lows += [mu - e for mu, e in zip(means, errs)]
        highs += [mu + e for mu, e in zip(means, errs)]
        ax.bar(x + (i - (len(trained) - 1) / 2) * w, means, w,
               color={'TransE': '#1f4e79', 'RotatE': '#6699cc'}[m], label=m,
               edgecolor='white', linewidth=0.4, zorder=2,
               yerr=errs, error_kw=dict(ecolor='#4a4a4a', elinewidth=0.8, capsize=2, capthick=0.8, zorder=3))
    ax.axhline(0, color='black', lw=0.8, zorder=1)
    ax.annotate('= Gemma', xy=(1.0, 0.0), xycoords=('axes fraction', 'data'),
                xytext=(4, 0), textcoords='offset points', ha='left', va='center',
                fontsize=6.5, color='#9a9a9a', annotation_clip=False)
    ymin = min(0, min(lows)); ymax = max(highs); pad = 0.06 * (ymax - ymin)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_ylabel('AUROC lift over Gemma prior (Δ)')
    ax.set_xticks(x)
    ax.set_xlim(-0.6, len(kgs) - 0.4)
    ax.tick_params(bottom=False)
    ax.set_xticklabels([KG_LABEL.get(k, k) for k in kgs])
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, ncol=len(trained), loc='lower center', bbox_to_anchor=(0.5, 1.01),
              handlelength=1.1, handletextpad=0.5, columnspacing=1.8, borderpad=0)
    ax.set_title('Lift of trained KG embeddings over the text-only name prior (Gemma)', loc='left', pad=20)
    fig.tight_layout()
    return fig


# ──────────────── 08 · heuristic vs embedding AUROC scatter (per strategy) ───────
def fig08_heuristic_scatter(single_csv, resampled_csv):
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Helvetica Neue", "Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7.5, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "text.color": "black", "axes.edgecolor": "black",
        "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
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
    pad = 0.03 * (vals.max() - vals.min())
    LO, HI = vals.min() - pad, vals.max() + pad

    fig, axes = plt.subplots(1, len(STRATEGIES), figsize=(7.4, 3.0), sharey=True)
    axes = np.atleast_1d(axes)
    for ax_i, strat in enumerate(STRATEGIES):
        ax = axes[ax_i]
        ax.set_axisbelow(True)
        ax.grid(True, color="#E8E8E8", lw=0.6, zorder=0)
        ax.plot([LO, HI], [LO, HI], ls=(0, (4, 3)), color="#BDBDBD", lw=0.8, zorder=1)
        for model_name, marker, ms in [("TransE", "o", 50), ("RotatE", "s", 38)]:
            sub = comp_df[(comp_df.model == model_name) & (comp_df.strategy == strat)].dropna(
                subset=["emb_auroc", "heuristic_auroc"])
            for _, r in sub.iterrows():
                ax.errorbar(r.heuristic_auroc, r.emb_auroc,
                            yerr=[[r.emb_auroc - r.emb_lo], [r.emb_hi - r.emb_auroc]],
                            fmt="none", ecolor=KG_COLORS.get(r.kg, "#888"),
                            elinewidth=0.8, capsize=2, capthick=0.8, alpha=0.7, zorder=4)
            ax.scatter(sub.heuristic_auroc, sub.emb_auroc,
                       c=[KG_COLORS.get(k, "#888") for k in sub.kg],
                       marker=marker, s=ms, zorder=5, edgecolors="white", linewidths=0.8)
        _sr = spearman_df[(spearman_df.model == "TransE") & (spearman_df.strategy == strat)]
        if len(_sr):
            r0 = _sr.iloc[0]
            ax.text(0.05, 0.95, f"TransE  ρ = {r0.spearman_rho:.2f}\np = {r0.p_value:.3f}",
                    transform=ax.transAxes, fontsize=7, va="top", color="black")
        ax.set_xlim(LO, HI); ax.set_ylim(LO, HI)
        ax.set_aspect("equal", adjustable="box")
        if ax_i == 0:
            ax.set_ylabel("Embedding AUROC")
        ax.set_title(strat, pad=6)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[-1].legend(
        handles=[Line2D([0], [0], marker="o", color="#555", lw=0, ms=6, label="TransE",
                        markeredgecolor="white", markeredgewidth=0.8),
                 Line2D([0], [0], marker="s", color="#555", lw=0, ms=5.5, label="RotatE",
                        markeredgecolor="white", markeredgewidth=0.8)],
        frameon=False, loc="lower right", handletextpad=0.4, borderpad=0.2)
    kg_handles = [Line2D([0], [0], marker="o", color=KG_COLORS[k], lw=0, ms=6,
                         markeredgecolor="white", markeredgewidth=0.8, label=KG_LABEL[k]) for k in KG_LABEL]
    fig.legend(handles=kg_handles, frameon=False, ncol=len(kg_handles), loc="upper center",
               bbox_to_anchor=(0.5, 1.0), handletextpad=0.3, columnspacing=1.1)
    fig.supxlabel("Heuristic AUROC", fontsize=8, y=0.02)
    fig.suptitle("Heuristic vs. embedding AUROC for drug–disease link prediction", y=1.10)
    fig.tight_layout()
    return fig
