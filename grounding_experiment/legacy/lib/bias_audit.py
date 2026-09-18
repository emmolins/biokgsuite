"""
Generic selection-bias audit: does a kept subset differ from the dropped remainder?

Dataset-agnostic. You supply a dataframe with a boolean `kept` column, a list of
categorical covariates and a list of continuous covariates; the functions run the
statistics and draw the figures. Nothing here is specific to any one benchmark.
"""
import numpy as np, pandas as pd
from scipy import stats as st
import matplotlib.pyplot as plt
try:                                    # shared house style (consistent across 09/09a/09b/09c)
    from . import figstyle
except ImportError:
    import figstyle
figstyle.apply()

KEPT_C, DROP_C = figstyle.PALETTE['covered'], figstyle.PALETTE['nokg']
_NEG, _NEUT = figstyle.PALETTE['nokg'], figstyle.PALETTE['pooled']


# ───────────────────────── statistics ─────────────────────────
def _odds_ratio(mask, kept):
    """OR of being kept given mask==True, with Haldane-corrected 95% CI + Fisher p."""
    a = (mask & kept).sum() + .5;  b = (mask & ~kept).sum() + .5
    c = (~mask & kept).sum() + .5; d = (~mask & ~kept).sum() + .5
    lor, se = np.log(a * d / (b * c)), np.sqrt(1/a + 1/b + 1/c + 1/d)
    p = st.fisher_exact(pd.crosstab(mask, kept))[1]
    return np.exp(lor), np.exp(lor - 1.96*se), np.exp(lor + 1.96*se), p


def run_bias_tests(df, kept_col, cat_vars, cont_vars, min_n=8):
    """Return a tidy results table: per categorical level (one-vs-rest OR) and per
    continuous covariate (Mann-Whitney). `min_n` skips tiny categorical levels."""
    kept = df[kept_col].astype(bool)
    rows = []
    for v in cat_vars:
        for lvl, n in df[v].value_counts().items():
            if n < min_n:
                continue
            orr, lo, hi, p = _odds_ratio(df[v] == lvl, kept)
            rows.append(dict(covariate=v, trait=f'{v} = {lvl}', kind='categorical',
                             n=int(n), effect=round(orr, 2), ci_low=round(lo, 2),
                             ci_high=round(hi, 2), p=round(p, 4)))
    for v in cont_vars:
        a, b = df.loc[kept, v].dropna(), df.loc[~kept, v].dropna()
        if len(a) and len(b):
            p = st.mannwhitneyu(a, b).pvalue
            rows.append(dict(covariate=v, trait=f'{v} (kept vs dropped median)', kind='continuous',
                             n=int(len(a) + len(b)), effect=f'{a.median():.1f} vs {b.median():.1f}',
                             ci_low=np.nan, ci_high=np.nan, p=round(p, 4)))
    return pd.DataFrame(rows)


# ───────────────────────── figures ─────────────────────────
def plot_forest(df, kept_col, traits, title='Odds of being kept, by trait'):
    """Forest plot of odds ratios. `traits` = {label: boolean Series}."""
    kept = df[kept_col].astype(bool)
    data = []
    for lbl, mask in traits.items():
        mask = mask.astype(bool)
        orr, lo, hi, p = _odds_ratio(mask, kept)
        data.append((lbl, orr, lo, hi, p))
    data.sort(key=lambda r: r[1])
    fig, ax = plt.subplots(figsize=(7.2, 0.55 * len(data) + 1.4))
    ax.axvline(1, color='#999', lw=1.1)
    for i, (lbl, orr, lo, hi, p) in enumerate(data):
        sig = p < 0.05
        c = _NEG if (sig and orr < 1) else (KEPT_C if sig else _NEUT)
        ax.plot([lo, hi], [i, i], color=c, lw=2.2, solid_capstyle='round', zorder=2)
        ax.scatter(orr, i, s=80, color=c, zorder=3, edgecolor='white', lw=1.1)
        ax.text(0.02, i, lbl, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=9)
        star = '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'n.s.'
        ax.text(hi * 1.05, i, star, va='center', fontsize=8.5, color=c, fontweight='bold' if sig else 'normal')
    ax.set_xscale('log'); ax.set_yticks([]); ax.set_ylim(-0.7, len(data) - 0.3)
    ax.set_xlabel('odds of being kept  (95% CI, log scale)')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    figstyle.title(ax, title, pad=12)
    fig.tight_layout(); return fig


def plot_distribution(df, kept_col, cont_var, title=None):
    """Strip + median comparison of a continuous covariate, kept vs dropped."""
    kept = df[kept_col].astype(bool)
    groups = [('kept', df[kept].copy(), KEPT_C), ('dropped', df[~kept].copy(), DROP_C)]
    fig, ax = plt.subplots(figsize=(7.2, 2.8)); rng = np.random.default_rng(3)
    for i, (lbl, sub, c) in enumerate(groups):
        v = sub[cont_var].dropna(); y = (1 - i) + rng.uniform(-0.16, 0.16, len(v))
        ax.scatter(v, y, s=26, color=c, alpha=0.55, edgecolor='white', lw=0.3)
        m = v.median(); ax.plot([m, m], [(1 - i) - 0.33, (1 - i) + 0.33], color=c, lw=3, solid_capstyle='round')
        ax.text(v.min(), 1 - i, f'  {lbl} (n={len(v)})', ha='left', va='center', fontsize=9, color=c, fontweight='bold')
    ax.set_yticks([]); ax.set_ylim(-0.7, 1.7)
    ax.set_xlabel(cont_var)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    figstyle.title(ax, title or f'{cont_var}: kept vs dropped', pad=10)
    fig.tight_layout(); return fig


def plot_composition(df, kept_col, cat_vars, title='Composition of kept vs dropped'):
    """Diverging paired bars (dropped ◀ | ▶ kept) of within-group share, one panel per
    categorical covariate. Panels share row spacing and are vertically centred."""
    kept = df[kept_col].astype(bool); grp = np.where(kept, 'kept', 'dropped')
    cts = []
    for v in cat_vars:
        ct = pd.crosstab(df[v], grp)
        ct = ct.loc[ct.sum(axis=1).sort_values().index]
        cts.append((v, ct))
    nmax = max(len(ct) for _, ct in cts)
    fig, axes = plt.subplots(1, len(cts), figsize=(3.3 * len(cts) + 0.6, 0.62 * nmax + 1.6))
    if len(cts) == 1: axes = [axes]
    n_k, n_d = int(kept.sum()), int((~kept).sum())
    for ax, (v, ct) in zip(axes, cts):
        k = ct['kept'] / ct['kept'].sum(); d = ct['dropped'] / ct['dropped'].sum()
        n = len(ct); off = (nmax - n) / 2.0; y = np.arange(n) + off
        ax.barh(y, -d.values, color=DROP_C, height=0.58); ax.barh(y, k.values, color=KEPT_C, height=0.58)
        ax.axvline(0, color='#bbb', lw=1.0)
        mx = max(k.max(), d.max())
        for yi, cat in zip(y, ct.index):
            ax.text(0, yi + 0.46, str(cat), ha='center', va='bottom', fontsize=9, color='#333')
            if d[cat] > 0.012: ax.text(-d[cat] - mx*0.03, yi, f'{d[cat]*100:.0f}%', ha='right', va='center', fontsize=7.8, color=DROP_C)
            if k[cat] > 0.012: ax.text(k[cat] + mx*0.03, yi, f'{k[cat]*100:.0f}%', ha='left', va='center', fontsize=7.8, color=KEPT_C)
        ax.set_xlim(-mx*1.32, mx*1.32); ax.set_ylim(-0.6, nmax - 0.3)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        ax.set_title(v, fontsize=10.5, fontweight='bold', pad=16)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=DROP_C, label=f'dropped (n={n_d})'), Patch(color=KEPT_C, label=f'kept (n={n_k})')],
               loc='upper center', ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 1.0))
    figstyle.suptitle(fig, title, y=1.08)
    fig.tight_layout(); return fig
