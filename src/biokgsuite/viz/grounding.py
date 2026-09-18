"""Figures for the grounding experiment, one function per display item.

Every function returns its Figure and saves nothing: the notebooks write
through ``biokgsuite.io.Outputs`` so each figure lands in the results manifest.
Graphs are coloured from ``KG_PALETTE``; models get the small fixed palette
below, assigned in the order the caller lists them. Display names come from
the ``*_label`` columns of the analysis tables or from a ``labels`` dict
(raw tag to display name), never from this module.
"""
from itertools import combinations

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from biokgsuite.grounding.analysis import (NONE_LABEL, cross_model_corr,
                                           mean_pairwise_corr, model_columns)
from biokgsuite.grounding.bias import odds_ratio
from biokgsuite.viz.style import (GRID_COLOR, KG_PALETTE, TEXT_COLOR, fig_size,
                                  panel_label, setup_style)

setup_style()

# Models are drawn beside graph colours in several panels, so their palette
# avoids the KG hues: no ocean blue, rust, sage, amber or plum.
MODEL_PALETTE = ['#3A4F7A', '#B07A2A', '#4E8A78', '#7A4E7A', '#8A6A4E', '#4E6A8A']

INK, MUTED, FAINT = '#1A1A1A', '#8C8C8C', '#BFBFBF'

# Pooled-vs-covered comparisons: covered takes the house blue because it is
# the arm that carries the claim, pooled is grey, and the no-KG baseline is
# the house rust so it reads as the reference line everywhere.
C_COVERED, C_POOLED, C_NOKG = '#2E6B8A', '#9A9A9A', '#C46A52'

# Bias audit: kept in the house blue, dropped in grey, a significant deficit
# in rust, a non-significant effect in neutral grey.
KEPT_COLOUR, DROPPED_COLOUR, NEG_COLOUR, NEUTRAL_COLOUR = '#2E6B8A', '#9A9A9A', '#C46A52', '#8C8C8C'

# Rank bins from best to worst: three blues fading to grey, then rust for a
# miss, so the eye reads "how much is dark" as "how often the drug was found".
RANK_RAMP = ['#2E6B8A', '#6D9BB5', '#A9C4D3', '#D9D9D9', '#C46A52']

# The evidence_agreement vocabulary of the ranking schema (prompt.RANKING_SCHEMA).
AGREEMENT_ORDER = ['consistent', 'conflicting', 'insufficient']
AGREEMENT_COLOURS = {'consistent': '#6A8A6B', 'conflicting': '#C49A3C',
                     'insufficient': '#C46A52'}
AGREEMENT_LABELS = {'consistent': 'Consistent\nevidence', 'conflicting': 'Conflicting',
                    'insufficient': 'Insufficient /\nmissing'}

# Ordinal ramp from light grey to the house blue, for "share of consistent
# evidence" and confidence colouring.
SEQ_CMAP = mcolors.LinearSegmentedColormap.from_list('grey_blue', [GRID_COLOR, '#2E6B8A'])


def _label(labels, key):
    if key == 'none':
        return (labels or {}).get('none', NONE_LABEL)
    return (labels or {}).get(key, key)


def _model_colours(models):
    return {m: MODEL_PALETTE[i % len(MODEL_PALETTE)] for i, m in enumerate(models)}


def _kg_colour(kg):
    return KG_PALETTE.get(kg, MUTED)


def _title(ax, title, subtitle=None, pad=10):
    ax.set_title(title, fontsize=11, fontweight='bold', loc='left', color=TEXT_COLOR, pad=pad)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=8.5, color='#777',
                va='bottom', ha='left')


def _panel(ax, letter, title, x=-0.13):
    panel_label(ax, letter, x=x, y=1.10, fontsize=9)
    ax.text(x + 0.02, 1.095, title, transform=ax.transAxes, fontsize=7.6,
            va='top', ha='left', color=INK)


def _bare(ax, spines=('top', 'right')):
    ax.grid(False)
    for side in spines:
        ax.spines[side].set_visible(False)


def _logx(ax, params):
    xs = sorted(set(float(p) for p in params))
    ax.set_xscale('log')
    ax.set_xticks(xs)
    ax.set_xticklabels([f'{p:g}B' for p in xs])
    ax.minorticks_off()
    ax.set_xlabel('Parameters (billions, log scale)')


def _first_seen(values):
    seen = []
    for v in values:
        if v not in seen:
            seen.append(v)
    return seen


def _pool_size(df):
    sizes = df['pool_size'].unique()
    if len(sizes) != 1:
        raise ValueError(f'expected one pool size, found {sorted(sizes)}')
    return int(sizes[0])


def fig_grounding_gain(arms_main, con_main, con_cov, model_order, kg_order, n_pairs, n_cov_pairs):
    """Two panels: the lift of every KG arm over no-KG per model, and every
    between-graph contrast on all pairs beside the same contrast on the pairs
    every graph covers.

    No gridlines, one reference rule per panel and direct labels instead of a
    legend. Panel b pairs each contrast on adjacent rows so the attenuation
    under coverage control reads as one movement rather than a scatter of
    intervals. The contrasts drawn are exactly those in ``con_main``; pass a
    filtered table to show fewer.
    """
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(fig_size(2)[0], 3.3),
                                   gridspec_kw={'width_ratios': [1.0, 1.3]})
    plt.subplots_adjust(left=0.085, right=0.985, bottom=0.155, top=0.86, wspace=0.50)

    width = 0.22
    offsets = {kg: (j - (len(kg_order) - 1) / 2) * width for j, kg in enumerate(kg_order)}
    lo_all, hi_all = [], []
    for kg in kg_order:
        for i, model in enumerate(model_order):
            hit = arms_main[(arms_main['model'] == model) & (arms_main['arm'] == kg)]
            if hit.empty:
                continue
            r = hit.iloc[0]
            x = i + offsets[kg]
            colour = _kg_colour(kg)
            axa.plot([x, x], [r['lift_lo'], r['lift_hi']], color=colour, lw=1.1,
                     solid_capstyle='round', zorder=2)
            axa.plot([x], [r['lift']], marker='o', ms=4.4, color=colour, mec='white',
                     mew=0.7, zorder=3)
            lo_all.append(r['lift_lo'])
            hi_all.append(r['lift_hi'])
            if i == 0:
                axa.annotate(r['arm_label'], (x, r['lift_hi']), xytext=(0, 5),
                             textcoords='offset points', ha='center', fontsize=6.3,
                             color=colour)
    model_labels = [arms_main.loc[arms_main['model'] == m, 'model_label'].iloc[0]
                    if (arms_main['model'] == m).any() else m for m in model_order]
    axa.set_xticks(range(len(model_order)))
    # Long model names collide at three per axis; break at the second hyphen
    # so each label sits on two lines under its own group.
    axa.set_xticklabels([lab.replace('-', '-\n', 1) if lab.count('-') >= 2 and len(lab) > 14 else lab
                         for lab in model_labels], fontsize=6.6)
    axa.set_ylabel('MRR gain over no-KG baseline', labelpad=3)
    if lo_all:
        span = max(hi_all) - min(lo_all)
        axa.set_ylim(min(lo_all) - 0.15 * span, max(hi_all) + 0.30 * span)
    axa.axhline(0, color=FAINT, lw=0.7, zorder=1)
    axa.set_xlim(-0.55, len(model_order) - 0.45)
    axa.tick_params(axis='x', length=0, pad=4)
    _bare(axa)
    _panel(axa, 'a', f'Grounding gain, all {n_pairs} pairs', x=-0.15)

    contrasts = _first_seen(con_main['contrast'])
    pair_off = 0.21
    y, ticks, tick_labels, heads = 0.0, [], [], []
    for model in model_order:
        head_label = (con_main.loc[con_main['model'] == model, 'model_label'].iloc[0]
                      if (con_main['model'] == model).any() else model)
        heads.append((y - 0.62, head_label))
        for con in contrasts:
            for frame, colour, off, lw, ms in ((con_main, FAINT, -pair_off, 1.0, 3.8),
                                               (con_cov, INK, pair_off, 1.3, 4.2)):
                hit = frame[(frame['model'] == model) & (frame['contrast'] == con)]
                if hit.empty:
                    continue
                r = hit.iloc[0]
                axb.plot([r['lo'], r['hi']], [y + off, y + off], color=colour, lw=lw,
                         solid_capstyle='round', zorder=2)
                axb.plot([r['diff']], [y + off], marker='o', ms=ms, zorder=3,
                         mfc=colour if (r['separates'] and colour == INK) else 'white',
                         mec=colour, mew=1.0)
            ticks.append(y)
            tick_labels.append(con.replace(' - ', ' − '))
            y += 1.30
        y += 0.95

    axb.axvline(0, color=INK, lw=0.8, zorder=1)
    axb.set_yticks(ticks)
    axb.set_yticklabels(tick_labels, fontsize=6.6)
    axb.set_ylim(y - 0.95 - 0.45, -1.25)
    axb.set_xlabel('Difference in MRR gain between graphs', labelpad=3)
    both = pd.concat([con_main, con_cov])
    if not both.empty:
        lo, hi = both['lo'].min(), both['hi'].max()
        span = hi - lo
        axb.set_xlim(min(lo, 0) - 0.15 * span, max(hi, 0) + 0.15 * span)
    axb.tick_params(axis='y', length=0, pad=2)
    _bare(axb, spines=('top', 'right', 'left'))
    x_head = axb.get_xlim()[0]
    for yh, label in heads:
        axb.text(x_head, yh, label, fontsize=6.9, fontweight='bold', color=INK,
                 va='center', ha='left')
    axb.annotate(f'pale, all {n_pairs} pairs   ·   dark, the {n_cov_pairs} covered by '
                 'every graph   ·   filled, excludes zero',
                 xy=(0.5, -0.155), xycoords='axes fraction', ha='center', va='top',
                 fontsize=6.2, color=MUTED)
    _panel(axb, 'b', 'Differences shrink under coverage control', x=-0.30)
    return fig


def fig_capacity(cap, chance_mrr, kg_order):
    """Model ladder: MRR per arm against parameter count, and the lift of each
    KG arm over no-KG. A flat no-KG line with a rising lift says the ability to
    use evidence, not prior knowledge, is what scales."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_size(2)[0], 2.75))
    plt.subplots_adjust(left=0.085, right=0.985, bottom=0.20, top=0.86, wspace=0.34)
    cx = cap.sort_values('params_b')
    xs = sorted(cx['params_b'].unique())

    for arm in list(kg_order) + ['none']:
        s = cx[cx['arm'] == arm].sort_values('params_b')
        if s.empty:
            continue
        ax1.plot(s['params_b'], s['mrr'], '-o', ms=3.6, lw=1.1, color=_kg_colour(arm),
                 label=s['arm_label'].iloc[0], zorder=3 if arm != 'none' else 2)
    ax1.axhline(chance_mrr, color=FAINT, ls=(0, (2, 2)), lw=0.7, zorder=1)
    ax1.text(xs[0], chance_mrr + 0.012, 'chance', fontsize=6, color=MUTED)
    _logx(ax1, xs)
    ax1.set_ylabel('Mean reciprocal rank', labelpad=3)
    ax1.grid(color=GRID_COLOR, lw=0.5)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, loc='upper left', handletextpad=0.4, labelspacing=0.3,
               borderaxespad=0.3)
    _panel(ax1, 'a', 'Ungrounded performance is flat')

    for arm in kg_order:
        s = cx[cx['arm'] == arm].sort_values('params_b')
        if s.empty:
            continue
        colour = _kg_colour(arm)
        ax2.errorbar(s['params_b'], s['lift'],
                     yerr=[s['lift'] - s['lift_lo'], s['lift_hi'] - s['lift']],
                     fmt='-o', ms=4.0, lw=1.1, color=colour, ecolor=colour,
                     elinewidth=0.9, capsize=0, zorder=3, label=s['arm_label'].iloc[0])
    ax2.axhline(0, color=INK, lw=0.7, zorder=1)
    _logx(ax2, xs)
    ax2.set_ylabel('MRR lift over no-KG', labelpad=3)
    ax2.grid(color=GRID_COLOR, lw=0.5)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, loc='upper left', handletextpad=0.4, labelspacing=0.3,
               borderaxespad=0.3)
    _panel(ax2, 'b', 'Ability to use evidence is what scales')
    return fig


def fig_lift_vs_size(scaling):
    """KG lift against model size, pooled over every KG arm and on covered
    arms only, with bootstrap intervals."""
    s = scaling.sort_values('params_b')
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.axhline(0, color=FAINT, ls='--', lw=0.8, zorder=1)
    for name, colour, marker, lo, hi, val in [
            ('covered arm', C_COVERED, 'D', s['lift_cov_lo'], s['lift_cov_hi'], s['lift_cov']),
            ('pooled (all KG arms)', C_POOLED, 'o', s['lift_pooled_lo'], s['lift_pooled_hi'],
             s['lift_pooled'])]:
        ax.vlines(s['params_b'], lo, hi, color=colour, lw=1, zorder=2)
        ax.plot(s['params_b'], val, marker=marker, ls='-', color=colour, lw=1.6, mfc='white',
                mec=colour, mew=1.4, ms=8, label=name, zorder=3)
    _logx(ax, s['params_b'])
    ax.set_ylabel('MRR lift (KG minus no-KG)')
    _title(ax, 'KG lift vs model size',
           'Mean reciprocal-rank lift over no-KG, pooled vs covered arm, 95% CI', pad=26)
    ax.legend(loc='lower right', handlelength=1.8)
    fig.tight_layout()
    return fig


def fig_mrr_vs_size(scaling):
    """Prior-only, pooled-KG and covered-arm MRR against model size. A flat
    no-KG line rules out memorisation as the source of the lift."""
    s = scaling.sort_values('params_b')
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for name, colour, marker, ls, val in [
            ('KG arm (covered)', C_COVERED, 'D', '-', s['mrr_kg_cov']),
            ('KG arm (pooled)', C_POOLED, 'o', '-', s['mrr_kg']),
            ('no-KG (prior only)', C_NOKG, 's', '--', s['mrr_none'])]:
        ax.plot(s['params_b'], val, marker=marker, ls=ls, color=colour, lw=1.6, mfc='white',
                mec=colour, mew=1.4, ms=7.5, label=name, zorder=3)
    _logx(ax, s['params_b'])
    ax.set_ylabel('Mean reciprocal rank')
    _title(ax, 'Prior-only vs KG-grounded MRR by model size',
           'no-KG, pooled KG and covered-arm MRR across the ladder', pad=26)
    ax.legend(loc='lower right', handlelength=2.2)
    fig.tight_layout()
    return fig


def fig_kg_by_model(df, params_b, kg_order, labels=None):
    """Grouped bars: each KG arm's MRR beside the no-KG baseline, one group per
    rung of the ladder, so a single graph driving the lift is visible."""
    models = [m for m in sorted(params_b, key=params_b.get) if (df['model'] == m).any()]
    arms = ['none'] + list(kg_order)
    mrr = df[df['model'].isin(models)].groupby(['model', 'arm'])['reciprocal_rank'].mean()
    x = np.arange(len(models))
    w = 0.8 / len(arms)
    fig, ax = plt.subplots(figsize=(max(6.0, 1.9 * len(models) + 2.0), 4.8))
    for i, arm in enumerate(arms):
        vals = [mrr.get((m, arm), np.nan) for m in models]
        ax.bar(x + (i - (len(arms) - 1) / 2) * w, vals, w, label=_label(labels, arm),
               color=_kg_colour(arm), edgecolor='white', linewidth=0.6, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([_label(labels, m) for m in models])
    ax.set_ylabel('Mean reciprocal rank')
    ax.set_ylim(0, None)
    _title(ax, 'MRR by knowledge graph and model size',
           'Mean reciprocal rank per KG arm vs the no-KG baseline', pad=26)
    ax.legend(ncol=len(arms), loc='upper left', columnspacing=1.3, handlelength=1.1)
    fig.tight_layout()
    return fig


def _rank_bins(ranks, pool_size):
    r = np.asarray(ranks)
    return [(r == 1).mean(), (r == 2).mean(), (r == 3).mean(),
            ((r >= 4) & (r <= pool_size)).mean(), (r > pool_size).mean()]


def fig_rank_distribution(df, models, arms, labels=None):
    """Stacked bars of where the true drug landed, one row per arm and one
    panel per model. The miss bin is rank > pool size."""
    pool_size = _pool_size(df)
    bins = ['Rank 1', 'Rank 2', 'Rank 3', f'Rank 4-{pool_size}', 'Missed']
    fig, axes = plt.subplots(1, len(models), figsize=(4.3 * len(models) + 0.2, 5.0),
                             sharex=True, squeeze=False)
    axes = axes[0]
    fig.subplots_adjust(left=0.05, right=0.985, top=0.80, bottom=0.11, wspace=0.08)
    for k, (ax, model) in enumerate(zip(axes, models)):
        for i, arm in enumerate(arms):
            sub = df[(df['model'] == model) & (df['arm'] == arm)]
            vals = _rank_bins(sub['rank'], pool_size) if len(sub) else [np.nan] * 5
            left, y = 0.0, len(arms) - 1 - i
            for v, colour in zip(vals, RANK_RAMP):
                if np.isnan(v):
                    continue
                ax.barh(y, v, left=left, height=0.72, color=colour, edgecolor='white', lw=0.6)
                if v > 0.05:
                    ax.text(left + v / 2, y, f'{v * 100:.0f}', ha='center', va='center',
                            fontsize=8.5,
                            color='white' if colour in (RANK_RAMP[0], RANK_RAMP[-1]) else '#333')
                left += v
        ax.set_yticks(range(len(arms)))
        ax.set_yticklabels([_label(labels, a) for a in arms][::-1] if k == 0 else [], fontsize=10)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xlabel('Fraction of queries', fontsize=10)
        ax.set_title(_label(labels, model), fontsize=12, fontweight='bold', pad=4)
        _bare(ax)
    handles = [Patch(color=c) for c in RANK_RAMP]
    fig.legend(handles, bins, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 0.96),
               fontsize=10, handlelength=1.1)
    return fig


def fig_evidence_quality(df):
    """Left: mean rank of the true drug by the model's own reading of the KG
    evidence. Right: hits@1 by self-reported confidence. KG arms only, because
    the no-KG arm has no evidence to agree with."""
    kg = df[df['arm'] != 'none']
    pool_size = _pool_size(df)
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(13.0, 5.4),
                                     gridspec_kw=dict(width_ratios=[1.15, 1]))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.84, bottom=0.13, wspace=0.55)

    for i, key in enumerate(AGREEMENT_ORDER):
        s = kg[kg['pos_agreement'] == key]
        y = len(AGREEMENT_ORDER) - 1 - i
        if len(s) == 0:
            continue
        colour = AGREEMENT_COLOURS[key]
        mean_rank = s['rank'].mean()
        top1 = (s['rank'] == 1).mean() * 100
        ax_l.plot([0, mean_rank], [y, y], color=colour, lw=4, solid_capstyle='round', zorder=2)
        ax_l.scatter(mean_rank, y, s=520, color=colour, zorder=3)
        ax_l.text(mean_rank, y, f'{mean_rank:.1f}', ha='center', va='center', color='white',
                  fontsize=10, fontweight='bold')
        ax_l.annotate(f'{top1:.0f}% ranked #1\n(n={len(s):,})', xy=(pool_size + 0.6, y),
                      ha='right', va='center', fontsize=10.5, color=colour, fontweight='bold')
    ax_l.set_yticks(range(len(AGREEMENT_ORDER)))
    ax_l.set_yticklabels([AGREEMENT_LABELS[k] for k in AGREEMENT_ORDER][::-1], fontsize=10.5)
    ax_l.set_xlim(0, pool_size + 3.5)
    ax_l.set_xticks(range(0, pool_size + 1))
    ax_l.set_ylim(-0.6, len(AGREEMENT_ORDER) - 0.4)
    ax_l.set_xlabel(f'Mean rank of the true drug (1 = best, {pool_size} = worst)', fontsize=10)
    ax_l.set_title('Mean rank of the true drug by KG evidence', fontsize=11.5,
                   fontweight='bold', loc='left', color=TEXT_COLOR, pad=6)
    _bare(ax_l)

    g = kg.groupby('pos_confidence')['hit1'].mean()
    x, yv = g.index.to_numpy(dtype=float), g.to_numpy(dtype=float)
    ax_r.grid(axis='y', color=GRID_COLOR, lw=0.8)
    ax_r.set_axisbelow(True)
    ax_r.plot(x, yv, color=MUTED, lw=1.6, zorder=2)
    ax_r.scatter(x, yv, s=300, c=[SEQ_CMAP(v) for v in yv], edgecolor='white', linewidth=1.3,
                 zorder=3)
    for xi, yi in zip(x, yv):
        ax_r.annotate(f'{yi * 100:.0f}%', (xi, yi), xytext=(0, 12), textcoords='offset points',
                      ha='center', fontsize=9.5, fontweight='bold', color='#333')
    ax_r.set_xticks([1, 2, 3, 4, 5])
    ax_r.set_xlim(0.6, 5.4)
    ax_r.set_ylim(-0.04, 1.1)
    ax_r.set_xlabel("Model's self-reported confidence (1-5)", fontsize=10)
    ax_r.set_ylabel('Hits@1', fontsize=10)
    ax_r.set_title('Hits@1 by model confidence', fontsize=11.5, fontweight='bold', loc='left',
                   color=TEXT_COLOR, pad=6)
    fig.suptitle('True-drug rank and hits@1 vs KG evidence', fontsize=14, fontweight='bold',
                 y=0.97, color=TEXT_COLOR)
    return fig


def fig_hardest_pairs(per, n=20, labels=None):
    """The ``n`` hardest pairs: pooled KG-arm MRR as a bar, per-model MRR as
    dots, so a pair that is hard for one model alone stands out."""
    cols = model_columns(per)
    colours = _model_colours(list(cols.values()))
    h = per.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.4, 0.34 * len(h) + 1.4))
    y = np.arange(len(h))
    ax.barh(y, h['mrr'], color=GRID_COLOR, height=0.64, zorder=1)
    for col, model in cols.items():
        ax.scatter(h[col], y, s=26, color=colours[model], zorder=3, edgecolor='white',
                   lw=0.4, label=_label(labels, model))
    ax.set_yticks(y)
    ax.set_yticklabels([f'{str(r.disease_name)[:40]}\n({str(r.drug_name)[:26]} · {r.area})'
                        for r in h.itertuples()], fontsize=6.8)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel('Mean MRR on KG arms (grey = pooled, dots = per model)')
    ax.legend(loc='lower right', fontsize=8, title='per-model MRR', title_fontsize=7.5)
    _title(ax, f'The {len(h)} hardest-to-rank drug-disease pairs', pad=10)
    _bare(ax)
    fig.tight_layout()
    return fig


def fig_difficulty_drivers(per):
    """Per-pair MRR against the model's mean self-reported confidence, coloured
    by the share of runs that read the KG evidence as consistent."""
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    sc = ax.scatter(per['pos_confidence'], per['mrr'], c=per['frac_consistent'], cmap=SEQ_CMAP,
                    s=46, edgecolor='white', lw=0.4, vmin=0, vmax=1)
    r = per[['pos_confidence', 'mrr']].corr().iloc[0, 1]
    ax.text(0.04, 0.96, f'Pearson r = {r:.2f}', transform=ax.transAxes, va='top', fontsize=10,
            color=INK, fontweight='bold')
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label('share of consistent KG evidence', fontsize=9)
    ax.set_xlabel("model's self-reported confidence (1-5)")
    ax.set_ylabel('mean MRR on KG arms')
    _title(ax, 'Per-pair MRR vs model confidence', pad=10)
    fig.tight_layout()
    return fig


def fig_difficulty_by_area(per, min_n=4):
    """Per-area difficulty: one dot per pair, a diamond at the area mean, areas
    sorted easiest to hardest from the top. Areas with fewer than ``min_n``
    pairs are left out because their mean is one or two pairs."""
    counts = per['area'].value_counts()
    keep = counts[counts >= min_n].index
    d = per[per['area'].isin(keep)]
    stat = d.groupby('area')['mrr'].agg(['mean', 'size']).sort_values('mean', ascending=False)
    order = stat.index.tolist()
    rng = np.random.default_rng(7)
    fig, ax = plt.subplots(figsize=(8.0, 0.62 * max(len(order), 1) + 1.3))
    overall = per['mrr'].mean()
    ax.axvline(overall, color=FAINT, ls='--', lw=1, zorder=1)
    for i, area in enumerate(order):
        y = len(order) - 1 - i
        sub = d[d['area'] == area]
        ax.scatter(sub['mrr'], y + rng.uniform(-0.15, 0.15, len(sub)), s=34, color='#A9C4D3',
                   alpha=0.9, edgecolor='white', lw=0.3, zorder=2)
        ax.scatter(sub['mrr'].mean(), y, s=150, color=C_COVERED, marker='D', zorder=4,
                   edgecolor='white', lw=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f'{a}  (n={int(stat.loc[a, "size"])})' for a in order[::-1]],
                       fontsize=9.5)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, max(len(order), 1) - 0.4)
    ax.set_xlabel('MRR on KG arms (small dot = pair, diamond = area mean)')
    ax.text(overall, len(order) - 0.5, 'overall mean', fontsize=8, color=MUTED, ha='center')
    _title(ax, 'Ranking difficulty by therapeutic area', pad=10)
    _bare(ax)
    fig.tight_layout()
    return fig


def fig_intrinsic_difficulty(per, labels=None):
    """Per-pair MRR of every model against every other, with the diagonal and
    Pearson r. Points hugging the diagonal mean difficulty is intrinsic to the
    pair rather than to the model."""
    cols = model_columns(per)
    pairs = list(combinations(list(cols), 2))
    fig, axes = plt.subplots(1, max(len(pairs), 1), figsize=(4.0 * max(len(pairs), 1), 4.2),
                             squeeze=False)
    axes = axes[0]
    for ax, (a, b) in zip(axes, pairs):
        d = per[[a, b]].dropna()
        ax.plot([0, 1], [0, 1], color=FAINT, ls='--', lw=1, zorder=1)
        ax.scatter(d[a], d[b], s=34, color=C_COVERED, alpha=0.7, edgecolor='white', lw=0.3,
                   zorder=2)
        r = d[a].corr(d[b]) if len(d) > 1 else np.nan
        ax.text(0.04, 0.93, f'r = {r:.2f}', transform=ax.transAxes, fontsize=11,
                fontweight='bold', color=INK)
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.set_aspect('equal')
        ax.set_xlabel(f'{_label(labels, cols[a])}, MRR')
        ax.set_ylabel(f'{_label(labels, cols[b])}, MRR')
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    r_mean = mean_pairwise_corr(cross_model_corr(per))
    fig.suptitle(f'Per-pair MRR agreement between models (mean pairwise r = {r_mean:.2f})',
                 fontsize=12, fontweight='bold', color=TEXT_COLOR, y=1.02)
    fig.tight_layout()
    return fig


def fig_coverage_fair(per):
    """MRR by the number of graphs covering the pair: pooled over every KG arm
    against covered arms only. A flat covered line under a rising pooled line
    means the coverage slope is the empty-arm effect, not pair difficulty."""
    d = per.dropna(subset=['n_kg_covered']).copy()
    d['n_kg_covered'] = d['n_kg_covered'].astype(int)
    d = d[d['n_kg_covered'] > 0]
    pooled = d.groupby('n_kg_covered')['mrr'].mean()
    covered = d.groupby('n_kg_covered')['mrr_covered'].mean()
    none = per['mrr_none'].mean()
    xs = sorted(d['n_kg_covered'].unique())
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.axhline(none, color=FAINT, ls=':', lw=1.2)
    if xs:
        ax.text(xs[-1] + 0.05, none, ' no-KG', va='center', fontsize=8.5, color=MUTED)
    ax.plot(xs, [pooled.get(k, np.nan) for k in xs], '-o', color=C_POOLED, lw=1.8, ms=9,
            mfc='white', mec=C_POOLED, mew=1.4, label='pooled (all KG arms)')
    ax.plot(xs, [covered.get(k, np.nan) for k in xs], '-D', color=C_COVERED, lw=1.8, ms=9,
            mfc='white', mec=C_COVERED, mew=1.4, label='covered arms only')
    for k in xs:
        if k in covered.index and not np.isnan(covered[k]):
            ax.text(k, covered[k] + 0.03, f'{covered[k]:.2f}', ha='center', fontsize=9,
                    color=C_COVERED, fontweight='bold')
        if k in pooled.index:
            ax.text(k, pooled[k] - 0.045, f'{pooled[k]:.2f}', ha='center', fontsize=9,
                    color=MUTED, fontweight='bold')
    ax.set_xticks(xs)
    ax.set_xticklabels([f'{k} KG' + ('s' if k > 1 else '') for k in xs])
    if xs:
        ax.set_xlim(xs[0] - 0.3, xs[-1] + 0.5)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel('number of KGs covering the pair')
    ax.set_ylabel('MRR on KG arms')
    ax.legend(loc='lower right')
    _title(ax, 'MRR by number of KGs covering the pair',
           'pooled over every arm vs covered arms only; the gap is the empty-arm effect',
           pad=26)
    fig.tight_layout()
    return fig


def _odds_rows(df, kept_col, traits):
    kept = df[kept_col].astype(bool)
    rows = [(label, *odds_ratio(mask.astype(bool), kept)) for label, mask in traits.items()]
    rows.sort(key=lambda r: r[1])
    return rows


def fig_bias_forest(df, kept_col, traits, title='Odds of being kept, by trait'):
    """Forest plot of the odds of being kept for every trait in ``traits``
    ({label: boolean Series}), sorted by effect, log scale."""
    data = _odds_rows(df, kept_col, traits)
    fig, ax = plt.subplots(figsize=(7.2, 0.55 * len(data) + 1.4))
    ax.axvline(1, color='#999', lw=1.1)
    for i, (label, orr, lo, hi, p) in enumerate(data):
        sig = p < 0.05
        colour = NEG_COLOUR if (sig and orr < 1) else (KEPT_COLOUR if sig else NEUTRAL_COLOUR)
        ax.plot([lo, hi], [i, i], color=colour, lw=2.2, solid_capstyle='round', zorder=2)
        ax.scatter(orr, i, s=80, color=colour, zorder=3, edgecolor='white', lw=1.1)
        ax.text(-0.02, i, label, transform=ax.get_yaxis_transform(), ha='right', va='center',
                fontsize=9)
        star = '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'n.s.'
        ax.text(hi * 1.05, i, star, va='center', fontsize=8.5, color=colour,
                fontweight='bold' if sig else 'normal')
    ax.set_xscale('log')
    # Plain numbers on the log axis: an odds ratio reads as 0.5 or 2, not 2x10^0.
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_yticks([])
    ax.set_ylim(-0.7, len(data) - 0.3)
    ax.set_xlabel('odds of being kept (95% CI, log scale)')
    _bare(ax)
    _title(ax, title, pad=12)
    fig.tight_layout()
    return fig


def fig_bias_distribution(df, kept_col, var, title=None):
    """Strip plot of one continuous covariate, kept above dropped, with a bar
    at each median."""
    kept = df[kept_col].astype(bool)
    groups = [('kept', df[kept], KEPT_COLOUR), ('dropped', df[~kept], DROPPED_COLOUR)]
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    rng = np.random.default_rng(3)
    for i, (label, sub, colour) in enumerate(groups):
        v = sub[var].dropna()
        y0 = 1 - i
        if len(v) == 0:
            continue
        ax.scatter(v, y0 + rng.uniform(-0.16, 0.16, len(v)), s=26, color=colour, alpha=0.55,
                   edgecolor='white', lw=0.3)
        med = v.median()
        ax.plot([med, med], [y0 - 0.33, y0 + 0.33], color=colour, lw=3, solid_capstyle='round')
        ax.text(v.min(), y0, f'  {label} (n={len(v)})', ha='left', va='center', fontsize=9,
                color=colour, fontweight='bold')
    ax.set_yticks([])
    ax.set_ylim(-0.7, 1.7)
    ax.set_xlabel(var)
    _bare(ax, spines=('top', 'right', 'left'))
    _title(ax, title or f'{var}: kept vs dropped', pad=10)
    fig.tight_layout()
    return fig


def fig_bias_composition(df, kept_col, cat_vars, title='Composition of kept vs dropped'):
    """Diverging paired bars (dropped to the left, kept to the right) of the
    within-group share of every level, one panel per categorical covariate."""
    kept = df[kept_col].astype(bool)
    grp = np.where(kept, 'kept', 'dropped')
    tables = []
    for var in cat_vars:
        ct = pd.crosstab(df[var], grp).reindex(columns=['kept', 'dropped'], fill_value=0)
        ct = ct.loc[ct.sum(axis=1).sort_values().index]
        tables.append((var, ct))
    n_max = max(len(ct) for _, ct in tables)
    fig, axes = plt.subplots(1, len(tables), figsize=(3.3 * len(tables) + 0.6, 0.62 * n_max + 1.6),
                             squeeze=False)
    axes = axes[0]
    n_k, n_d = int(kept.sum()), int((~kept).sum())
    for ax, (var, ct) in zip(axes, tables):
        k = ct['kept'] / max(ct['kept'].sum(), 1)
        d = ct['dropped'] / max(ct['dropped'].sum(), 1)
        n = len(ct)
        y = np.arange(n) + (n_max - n) / 2.0
        ax.barh(y, -d.to_numpy(), color=DROPPED_COLOUR, height=0.58)
        ax.barh(y, k.to_numpy(), color=KEPT_COLOUR, height=0.58)
        ax.axvline(0, color='#bbb', lw=1.0)
        mx = max(k.max(), d.max(), 1e-9)
        for yi, cat in zip(y, ct.index):
            ax.text(0, yi + 0.46, str(cat), ha='center', va='bottom', fontsize=9, color='#333')
            if d[cat] > 0.012:
                ax.text(-d[cat] - mx * 0.03, yi, f'{d[cat] * 100:.0f}%', ha='right',
                        va='center', fontsize=7.8, color=DROPPED_COLOUR)
            if k[cat] > 0.012:
                ax.text(k[cat] + mx * 0.03, yi, f'{k[cat] * 100:.0f}%', ha='left', va='center',
                        fontsize=7.8, color=KEPT_COLOUR)
        ax.set_xlim(-mx * 1.32, mx * 1.32)
        ax.set_ylim(-0.6, n_max - 0.3)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(var, fontsize=10.5, fontweight='bold', pad=16)
    fig.legend(handles=[Patch(color=DROPPED_COLOUR, label=f'dropped (n={n_d})'),
                        Patch(color=KEPT_COLOUR, label=f'kept (n={n_k})')],
               loc='upper center', ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(title, fontsize=12, fontweight='bold', color=TEXT_COLOR, y=1.08)
    fig.tight_layout()
    return fig
