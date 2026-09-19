"""Figures for the grounding experiment, one function per display item.

Every function returns its Figure and saves nothing: the notebooks write
through ``biokgsuite.io.Outputs`` so each figure lands in the results manifest.
Graphs are coloured from ``KG_PALETTE``; models get the small fixed palette
below, assigned in the order the caller lists them, doubled by a marker shape
so a figure that compares models still reads in greyscale. The no-KG arm and
every other baseline is ``NEUTRAL`` grey. Display names come from the
``*_label`` columns of the analysis tables or from a ``labels`` dict (raw tag
to display name), never from this module.
"""
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from biokgsuite.grounding.analysis import (NONE_LABEL, cross_model_corr,
                                           mean_pairwise_corr, model_columns)
from biokgsuite.grounding.bias import odds_ratio
from biokgsuite.viz.primitives import (facet_grid, human_label, panel_letters,
                                       ref_line, truncate)
from biokgsuite.viz.style import (DOUBLE_COL_W, FAINT, HEATMAP_CMAP, KG_PALETTE, MUTED, NEUTRAL, SINGLE_COL_W, TEXT_COLOR, clean_ax,
                                  facet_title, fig_size, panel_label, setup_style)

setup_style()

INK = TEXT_COLOR

# Models are not graphs, so they never borrow a KG hue. Three muted hues
# (periwinkle, ochre, wine) that pass the dataviz palette validator among
# themselves under protan, deutan and normal vision, and sit as far from the
# Okabe-Ito graph hues as the remaining hue circle allows. A marker shape is
# always drawn with the colour, so the assignment survives greyscale print.
MODEL_PALETTE = ['#8085C8', '#977D30', '#84446A']
MODEL_MARKERS = ['o', 's', '^']

# Rank bins from best to worst: three steps of the magnitude ramp, dark to
# light, then faint grey for the tail and a white cell with a hairline edge
# for a miss, so "how much is dark" reads as "how often the drug was found"
# without a warning colour for what is only the last bin.
RANK_RAMP = [HEATMAP_CMAP(0.92), HEATMAP_CMAP(0.62), HEATMAP_CMAP(0.36), FAINT, 'white']

# The evidence_agreement vocabulary of the ranking schema (prompt.RANKING_SCHEMA).
AGREEMENT_ORDER = ['consistent', 'conflicting', 'insufficient']
AGREEMENT_LABELS = {'consistent': 'Consistent', 'conflicting': 'Conflicting',
                    'insufficient': 'Insufficient or missing'}

BAR_W = 0.6


def _label(labels, key):
    if key == 'none':
        return (labels or {}).get('none', NONE_LABEL)
    return (labels or {}).get(key, key)


def _model_colours(models):
    return {m: MODEL_PALETTE[i % len(MODEL_PALETTE)] for i, m in enumerate(models)}


def _model_markers(models):
    return {m: MODEL_MARKERS[i % len(MODEL_MARKERS)] for i, m in enumerate(models)}


def _kg_colour(kg):
    return KG_PALETTE.get(kg, NEUTRAL)


def _two_lines(label, limit=12):
    """Break a long model name at its second hyphen so three names fit under
    one axis without rotating; short names stay on one line."""
    lab = str(label)
    if len(lab) <= limit or lab.count('-') < 2:
        return lab
    first = lab.find('-', lab.find('-') + 1)
    return lab[:first] + '-\n' + lab[first + 1:]


def _vref(ax, x, label, pad_pt=2):
    """Vertical counterpart of ``ref_line``: grey dashed with a 6 pt muted
    label just above the top of the axes."""
    ax.axvline(x, color=FAINT, ls='--', lw=0.6, zorder=1)
    ax.annotate(label, xy=(x, 1), xycoords=('data', 'axes fraction'),
                xytext=(0, pad_pt), textcoords='offset points',
                ha='center', va='bottom', fontsize=6, color=MUTED, annotation_clip=False)


def _logx(ax, params):
    """Parameter count on a log axis with the rungs themselves as ticks."""
    xs = sorted(set(float(p) for p in params))
    ax.set_xscale('log')
    ax.set_xticks(xs)
    ax.set_xticklabels([f'{p:g}B' for p in xs])
    ax.minorticks_off()
    ax.set_xlabel('Parameters (billions)')


def _end_labels_or_legend(ax, ends, loc='best', min_gap=0.07):
    """Direct labels at the line ends when they are a text height apart,
    otherwise a 6 pt frameless legend. ``ends`` is {label: (x, y, colour)}."""
    ys = sorted(v[1] for v in ends.values() if v[1] == v[1])
    lo, hi = ax.get_ylim()
    gaps = np.diff(ys) if len(ys) > 1 else np.array([np.inf])
    if ends and gaps.min() >= min_gap * (hi - lo):
        x0, x1 = ax.get_xlim()
        if ax.get_xscale() == 'log':
            ax.set_xlim(x0, x1 * (x1 / x0) ** 0.22)
        else:
            ax.set_xlim(x0, x1 + 0.22 * (x1 - x0))
        for label, (x, y, colour) in ends.items():
            if y == y:
                ax.annotate(label, (x, y), xytext=(4, 0), textcoords='offset points',
                            ha='left', va='center', fontsize=6, color=colour)
    elif len(ends) > 1:
        ax.legend(loc=loc)


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


def _dark(colour):
    r, g, b = plt.matplotlib.colors.to_rgb(colour)
    return 0.299 * r + 0.587 * g + 0.114 * b < 0.5


def fig_grounding_gain(arms_main, con_main, con_cov, model_order, kg_order, n_pairs, n_cov_pairs):
    """Two panels: the lift of every KG arm over no-KG per model, and every
    between-graph contrast on all pairs beside the same contrast on the pairs
    every graph can evaluate (holds and has evidence for).

    Panel b pairs each contrast on adjacent rows so the attenuation under
    evidence control reads as one movement rather than a scatter of
    intervals. The contrasts drawn are exactly those in ``con_main``; pass a
    filtered table to show fewer.
    """
    contrasts = _first_seen(con_main['contrast'])
    n_rows = len(model_order) * (len(contrasts) + 1)
    # Nine contrast rows (three models, three graphs) fit the standard row
    # height; fifteen contrasts per model need the panel to grow with them.
    height = max(fig_size(2)[1] * 1.3, 0.13 * n_rows + 0.7)
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(fig_size(2)[0], height),
                                   gridspec_kw={'width_ratios': [1.0, 1.3]})
    plt.subplots_adjust(left=0.08, right=0.985, bottom=0.6 / height, top=1 - 0.22 / height,
                        wspace=0.5)

    width = min(0.22, 0.8 / max(len(kg_order), 1))
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
            axa.plot([x, x], [r['lift_lo'], r['lift_hi']], color=colour, lw=0.75,
                     solid_capstyle='butt', zorder=2)
            axa.plot([x], [r['lift']], marker='o', ms=3.5, color=colour, mec='white',
                     mew=0.5, zorder=3)
            lo_all.append(r['lift_lo'])
            hi_all.append(r['lift_hi'])
            if i == 0:
                # The graph is named once, above its first interval, in its own
                # colour: a direct label in place of a legend.
                axa.annotate(r['arm_label'], (x, r['lift_hi']), xytext=(0, 3),
                             textcoords='offset points', ha='center', va='bottom',
                             fontsize=6, color=colour)
    model_labels = [arms_main.loc[arms_main['model'] == m, 'model_label'].iloc[0]
                    if (arms_main['model'] == m).any() else m for m in model_order]
    axa.set_xticks(range(len(model_order)))
    axa.set_xticklabels([_two_lines(lab) for lab in model_labels])
    if lo_all:
        # The axis always includes zero: a gain is read against no gain.
        span = max(hi_all) - min(lo_all + [0])
        axa.set_ylim(min(lo_all + [0]) - 0.05 * span, max(hi_all) + 0.30 * span)
    axa.set_xlim(-0.55, len(model_order) - 0.45)
    ref_line(axa, 0, 'no gain')
    clean_ax(axa, ylabel='MRR gain over no-KG')
    axa.tick_params(axis='x', length=0)

    pair_off = 0.21
    y, ticks, tick_labels, heads = 0.0, [], [], []
    for model in model_order:
        head_label = (con_main.loc[con_main['model'] == model, 'model_label'].iloc[0]
                      if (con_main['model'] == model).any() else model)
        heads.append((y - 0.62, head_label))
        for con in contrasts:
            for frame, colour, off, ms in ((con_main, FAINT, -pair_off, 3.0),
                                           (con_cov, INK, pair_off, 3.5)):
                hit = frame[(frame['model'] == model) & (frame['contrast'] == con)]
                if hit.empty:
                    continue
                r = hit.iloc[0]
                axb.plot([r['lo'], r['hi']], [y + off, y + off], color=colour, lw=0.75,
                         solid_capstyle='butt', zorder=2)
                axb.plot([r['diff']], [y + off], marker='o', ms=ms, zorder=3,
                         mfc=colour if (r['separates'] and colour == INK) else 'white',
                         mec=colour, mew=0.75)
            ticks.append(y)
            tick_labels.append(con.replace(' - ', ' − '))
            y += 1.30
        y += 0.95

    axb.axvline(0, color=FAINT, ls='--', lw=0.6, zorder=1)
    axb.set_yticks(ticks)
    axb.set_yticklabels(tick_labels)
    axb.set_ylim(y - 0.95 - 0.45, -1.25)
    both = pd.concat([con_main, con_cov])
    if not both.empty:
        lo, hi = both['lo'].min(), both['hi'].max()
        span = hi - lo
        axb.set_xlim(min(lo, 0) - 0.15 * span, max(hi, 0) + 0.15 * span)
    clean_ax(axb, xlabel='Difference in MRR gain between graphs')
    axb.spines['left'].set_visible(False)
    axb.tick_params(axis='y', length=0, pad=2)
    # Model headings sit at the left edge of the panel at the same size as the
    # tick labels, regular weight: only panel letters are bold.
    x_head = axb.get_xlim()[0]
    for yh, label in heads:
        axb.text(x_head, yh, label, fontsize=6.5, color=INK, va='center', ha='left')
    axb.annotate(f'pale, all {n_pairs} pairs; dark, the {n_cov_pairs} evaluable in every graph;\n'
                 'filled, interval excludes zero',
                 xy=(0, 0), xycoords='axes fraction', xytext=(0, -26),
                 textcoords='offset points', ha='left', va='top', fontsize=6,
                 color=MUTED, annotation_clip=False)
    panel_label(axa, 'a', x=-0.28, y=1.02)
    panel_label(axb, 'b', x=-0.30, y=1.02)
    return fig


def fig_capacity(cap, chance_mrr, kg_order):
    """Model ladder: MRR per arm against parameter count, and the lift of each
    KG arm over no-KG. A flat no-KG line with a rising lift says the ability to
    use evidence, not prior knowledge, is what scales."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=fig_size(2, 1), constrained_layout=True)
    cx = cap.sort_values('params_b')
    xs = sorted(cx['params_b'].unique())

    for arm in list(kg_order) + ['none']:
        s = cx[cx['arm'] == arm].sort_values('params_b')
        if s.empty:
            continue
        none = arm == 'none'
        ax1.plot(s['params_b'], s['mrr'], ls='--' if none else '-', marker='s' if none else 'o',
                 ms=3.5, lw=1, color=_kg_colour(arm), label=s['arm_label'].iloc[0],
                 zorder=2 if none else 3)
    ref_line(ax1, chance_mrr, 'chance')
    _logx(ax1, xs)
    clean_ax(ax1, ylabel='Mean reciprocal rank')
    ax1.legend(loc='best', ncol=1 if len(kg_order) <= 3 else 2)

    for arm in kg_order:
        s = cx[cx['arm'] == arm].sort_values('params_b')
        if s.empty:
            continue
        colour = _kg_colour(arm)
        ax2.errorbar(s['params_b'], s['lift'],
                     yerr=[s['lift'] - s['lift_lo'], s['lift_hi'] - s['lift']],
                     fmt='-o', ms=3.5, lw=1, color=colour, ecolor=colour,
                     elinewidth=0.75, capsize=0, zorder=3, label=s['arm_label'].iloc[0])
    ref_line(ax2, 0, 'no gain')
    _logx(ax2, xs)
    clean_ax(ax2, ylabel='MRR gain over no-KG')
    ax2.legend(loc='best', ncol=1 if len(kg_order) <= 3 else 2)
    panel_letters(fig, [ax1, ax2])
    return fig


def fig_lift_vs_size(scaling):
    """KG lift against model size, pooled over every KG arm and on grounded
    arms only (the graph had evidence for the positive), with bootstrap
    intervals. Neither series is a graph, so the grounded arm is ink and the
    pooled arm grey."""
    s = scaling.sort_values('params_b')
    fig, ax = plt.subplots(figsize=fig_size(1, 1), constrained_layout=True)
    ends = {}
    for name, colour, marker, lo, hi, val in [
            ('Grounded arms', INK, 'D', s['lift_grounded_lo'], s['lift_grounded_hi'], s['lift_grounded']),
            ('Pooled, all KG arms', NEUTRAL, 'o', s['lift_pooled_lo'], s['lift_pooled_hi'],
             s['lift_pooled'])]:
        ax.vlines(s['params_b'], lo, hi, color=colour, lw=0.75, zorder=2)
        ax.plot(s['params_b'], val, marker=marker, ls='-', color=colour, lw=1, mfc='white',
                mec=colour, mew=0.75, ms=3.5, label=name, zorder=3)
        if len(s):
            ends[name] = (s['params_b'].iloc[-1], val.iloc[-1], colour)
    ref_line(ax, 0, 'no gain')
    _logx(ax, s['params_b'])
    clean_ax(ax, ylabel='MRR gain over no-KG')
    _end_labels_or_legend(ax, ends)
    return fig


def fig_mrr_vs_size(scaling):
    """No-KG, pooled-KG and grounded-arm MRR against model size. A flat
    no-KG line rules out memorisation as the source of the lift."""
    s = scaling.sort_values('params_b')
    fig, ax = plt.subplots(figsize=fig_size(1, 1), constrained_layout=True)
    ends = {}
    for name, colour, marker, ls, val in [
            ('Grounded arms', INK, 'D', '-', s['mrr_kg_grounded']),
            ('Pooled, all KG arms', INK, 'o', ':', s['mrr_kg']),
            ('No KG', NEUTRAL, 's', '--', s['mrr_none'])]:
        ax.plot(s['params_b'], val, marker=marker, ls=ls, color=colour, lw=1, mfc='white',
                mec=colour, mew=0.75, ms=3.5, label=name, zorder=3)
        if len(s):
            ends[name] = (s['params_b'].iloc[-1], val.iloc[-1], colour)
    _logx(ax, s['params_b'])
    clean_ax(ax, ylabel='Mean reciprocal rank')
    _end_labels_or_legend(ax, ends)
    return fig


def fig_kg_by_model(df, params_b, kg_order, labels=None):
    """Grouped bars: each KG arm's MRR beside the no-KG baseline, one group per
    rung of the ladder, so a single graph driving the lift is visible."""
    models = [m for m in sorted(params_b, key=params_b.get) if (df['model'] == m).any()]
    arms = ['none'] + list(kg_order)
    mrr = df[df['model'].isin(models)].groupby(['model', 'arm'])['reciprocal_rank'].mean()
    x = np.arange(len(models))
    w = BAR_W / len(arms)
    cols = 1 if len(models) * len(arms) <= 16 else 2
    fig, ax = plt.subplots(figsize=fig_size(cols, 1), constrained_layout=True)
    top = 0
    for i, arm in enumerate(arms):
        vals = [mrr.get((m, arm), np.nan) for m in models]
        top = max(top, np.nanmax(vals) if np.isfinite(vals).any() else 0)
        ax.bar(x + (i - (len(arms) - 1) / 2) * w, vals, w, label=_label(labels, arm),
               color=_kg_colour(arm), zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([_two_lines(_label(labels, m)) for m in models])
    ax.set_xlim(-0.6, len(models) - 0.4)
    # Head room for a frameless legend above the tallest bar.
    ax.set_ylim(0, top * 1.3 if top else 1)
    clean_ax(ax, ylabel='Mean reciprocal rank')
    ax.tick_params(axis='x', length=0)
    ax.legend(loc='upper left', ncol=min(len(arms), 4))
    return fig


def _rank_bins(ranks, pool_size):
    r = np.asarray(ranks)
    return [(r == 1).mean(), (r == 2).mean(), (r == 3).mean(),
            ((r >= 4) & (r <= pool_size)).mean(), (r > pool_size).mean()]


def fig_rank_distribution(df, models, arms, labels=None):
    """Stacked bars of where the true drug landed, one row per arm and one
    panel per model. The miss bin is rank > pool size."""
    pool_size = _pool_size(df)
    bins = ['Rank 1', 'Rank 2', 'Rank 3', f'Rank 4–{pool_size}', 'Missed']
    row_h = max(fig_size(1)[1], 0.16 * len(arms) + 0.7)
    fig, axes, ncol = facet_grid(len(models), row_h=row_h)
    for k, (ax, model) in enumerate(zip(axes, models)):
        for i, arm in enumerate(arms):
            sub = df[(df['model'] == model) & (df['arm'] == arm)]
            vals = _rank_bins(sub['rank'], pool_size) if len(sub) else [np.nan] * 5
            left, y = 0.0, len(arms) - 1 - i
            for v, colour in zip(vals, RANK_RAMP):
                if np.isnan(v):
                    continue
                miss = colour == 'white'
                ax.barh(y, v, left=left, height=BAR_W, color=colour,
                        edgecolor=NEUTRAL if miss else 'white', lw=0.4, zorder=2)
                # A two-digit count needs about 8 pt; a segment narrower than
                # 7% of a 55 mm panel cannot hold it.
                if v > 0.07:
                    ax.text(left + v / 2, y, f'{v * 100:.0f}', ha='center', va='center',
                            fontsize=6, color='white' if _dark(colour) else TEXT_COLOR)
                left += v
        ax.set_yticks(range(len(arms)))
        ax.set_yticklabels([_label(labels, a) for a in arms][::-1] if k % ncol == 0 else [])
        ax.set_ylim(-0.6, len(arms) - 0.4)
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1])
        clean_ax(ax, xlabel='Fraction of queries' if k >= len(models) - ncol else '')
        ax.tick_params(axis='y', length=0)
        facet_title(ax, _label(labels, model))
    handles = [Patch(facecolor=c, edgecolor=NEUTRAL if c == 'white' else 'none', lw=0.4)
               for c in RANK_RAMP]
    fig.legend(handles, bins, ncol=5, loc='outside lower center')
    if len(models) > 1:
        panel_letters(fig, axes)
    return fig


def fig_evidence_quality(df):
    """Left: mean rank of the true drug by the model's own reading of the KG
    evidence. Right: hits@1 by self-reported confidence. KG arms only, because
    the no-KG arm has no evidence to agree with."""
    kg = df[df['arm'] != 'none']
    pool_size = _pool_size(df)
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=fig_size(2, 1), constrained_layout=True,
                                     gridspec_kw=dict(width_ratios=[1.15, 1]))

    for i, key in enumerate(AGREEMENT_ORDER):
        s = kg[kg['pos_agreement'] == key]
        y = len(AGREEMENT_ORDER) - 1 - i
        if len(s) == 0:
            continue
        mean_rank = s['rank'].mean()
        top1 = (s['rank'] == 1).mean() * 100
        ax_l.plot([0, mean_rank], [y, y], color=INK, lw=1, solid_capstyle='butt', zorder=2)
        ax_l.plot([mean_rank], [y], marker='o', ms=3.5, color=INK, zorder=3)
        ax_l.annotate(f'{mean_rank:.1f}', (mean_rank, y), xytext=(0, 4),
                      textcoords='offset points', ha='center', va='bottom', fontsize=6,
                      color=TEXT_COLOR)
        ax_l.annotate(f'{top1:.0f}% ranked first (n = {len(s):,})', xy=(1, y),
                      xycoords=('axes fraction', 'data'), xytext=(-2, 0),
                      textcoords='offset points', ha='right', va='center', fontsize=6,
                      color=MUTED)
    ax_l.set_yticks(range(len(AGREEMENT_ORDER)))
    ax_l.set_yticklabels([AGREEMENT_LABELS[k] for k in AGREEMENT_ORDER][::-1])
    ax_l.set_xlim(0, pool_size + 3.5)
    ax_l.set_xticks(range(0, pool_size + 1))
    ax_l.set_ylim(-0.6, len(AGREEMENT_ORDER) - 0.4)
    clean_ax(ax_l, xlabel=f'Mean rank of the true drug (1 best, {pool_size} worst)')
    ax_l.tick_params(axis='y', length=0)
    facet_title(ax_l, 'By reported evidence agreement')

    g = kg.groupby('pos_confidence')['hit1'].mean()
    x, yv = g.index.to_numpy(dtype=float), g.to_numpy(dtype=float)
    ax_r.plot(x, yv, color=INK, lw=1, marker='o', ms=3.5, zorder=3)
    ax_r.set_xticks([1, 2, 3, 4, 5])
    ax_r.set_xlim(0.6, 5.4)
    ax_r.set_ylim(0, 1.0)
    clean_ax(ax_r, xlabel='Self-reported confidence (1 low, 5 high)', ylabel='Hits@1')
    facet_title(ax_r, 'By reported confidence')
    panel_letters(fig, [ax_l, ax_r])
    return fig


def fig_hardest_pairs(per, n=20, labels=None):
    """The ``n`` hardest pairs: pooled KG-arm MRR as a bar, per-model MRR as
    dots, so a pair that is hard for one model alone stands out."""
    cols = model_columns(per)
    models = list(cols.values())
    colours, markers = _model_colours(models), _model_markers(models)
    h = per.head(n).iloc[::-1]
    # Forty-character pair names take a third of the double column at 6 pt.
    fig, ax = plt.subplots(figsize=(DOUBLE_COL_W, 0.13 * len(h) + 0.55), constrained_layout=True)
    y = np.arange(len(h))
    ax.barh(y, h['mrr'], color=FAINT, height=BAR_W, zorder=1)
    for col, model in cols.items():
        ax.scatter(h[col], y, s=12, color=colours[model], marker=markers[model], zorder=3,
                   edgecolor='white', lw=0.3, label=_label(labels, model))
    ax.set_yticks(y)
    ax.set_yticklabels([truncate(f'{r.disease_name} ({r.drug_name})', 40) for r in h.itertuples()])
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.6, len(h) - 0.4)
    clean_ax(ax, xlabel='Mean reciprocal rank on KG arms (bar, pooled; dots, per model)')
    ax.tick_params(axis='y', length=0)
    ax.legend(loc='lower right')
    return fig


def fig_difficulty_drivers(per):
    """Per-pair MRR against the model's mean self-reported confidence, coloured
    by the share of runs that read the KG evidence as consistent."""
    fig, ax = plt.subplots(figsize=(SINGLE_COL_W, 2.4), constrained_layout=True)
    sc = ax.scatter(per['pos_confidence'], per['mrr'], c=per['frac_consistent'],
                    cmap=HEATMAP_CMAP, s=10, alpha=0.8, lw=0, vmin=0, vmax=1, zorder=2)
    r = per[['pos_confidence', 'mrr']].corr().iloc[0, 1]
    ax.text(0.03, 0.97, f'Pearson r = {r:.2f}', transform=ax.transAxes, va='top', ha='left',
            fontsize=6, color=MUTED)
    cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cb.ax.tick_params(labelsize=6, length=2, width=0.5, colors=TEXT_COLOR)
    cb.outline.set_linewidth(0.5)
    cb.set_label('Share of consistent evidence', fontsize=6, color=TEXT_COLOR, labelpad=3)
    # Confidence is reported on a fixed 1 to 5 scale, so the axis shows it whole.
    ax.set_xlim(0.8, 5.2)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_ylim(0, 1.02)
    clean_ax(ax, xlabel='Self-reported confidence (1 low, 5 high)', ylabel='Mean MRR on KG arms')
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
    fig, ax = plt.subplots(figsize=(SINGLE_COL_W, 0.16 * max(len(order), 1) + 0.7),
                           constrained_layout=True)
    overall = per['mrr'].mean()
    for i, area in enumerate(order):
        y = len(order) - 1 - i
        sub = d[d['area'] == area]
        ax.scatter(sub['mrr'], y + rng.uniform(-0.15, 0.15, len(sub)), s=8, color=NEUTRAL,
                   alpha=0.6, lw=0, zorder=2, label='Pair' if i == 0 else None)
        ax.scatter(sub['mrr'].mean(), y, s=18, color=INK, marker='D', zorder=4,
                   edgecolor='white', lw=0.4, label='Area mean' if i == 0 else None)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f'{a} (n = {int(stat.loc[a, "size"])})' for a in order[::-1]])
    ax.set_xlim(0, 1.02)
    ax.set_ylim(-0.6, max(len(order), 1) - 0.4)
    _vref(ax, overall, 'overall mean')
    clean_ax(ax, xlabel='Mean reciprocal rank on KG arms')
    ax.tick_params(axis='y', length=0)
    if order:
        ax.legend(loc='lower right')
    return fig


def fig_intrinsic_difficulty(per, labels=None):
    """Per-pair MRR of every model against every other, with the diagonal and
    Pearson r. Points hugging the diagonal mean difficulty is intrinsic to the
    pair rather than to the model. The mean pairwise r is computed for the
    caption and the notebook, not drawn."""
    cols = model_columns(per)
    pairs = list(combinations(list(cols), 2))
    n = max(len(pairs), 1)
    # Square panels: three of them across the double column are 55 mm wide.
    row_h = (SINGLE_COL_W if n == 1 else DOUBLE_COL_W / min(n, 3)) + 0.25
    fig, axes, ncol = facet_grid(n, row_h=row_h)
    for ax, (a, b) in zip(axes, pairs):
        d = per[[a, b]].dropna()
        ax.plot([0, 1], [0, 1], color=FAINT, ls='--', lw=0.6, zorder=1)
        ax.scatter(d[a], d[b], s=8, color=INK, alpha=0.6, lw=0, zorder=2)
        r = d[a].corr(d[b]) if len(d) > 1 else np.nan
        ax.text(0.04, 0.96, f'r = {r:.2f}', transform=ax.transAxes, fontsize=6, va='top',
                color=MUTED)
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)
        ax.set_aspect('equal')
        clean_ax(ax, xlabel=f'{_label(labels, cols[a])}, MRR',
                 ylabel=f'{_label(labels, cols[b])}, MRR')
    mean_pairwise_corr(cross_model_corr(per))
    if len(pairs) > 1:
        panel_letters(fig, axes)
    return fig


def fig_coverage_fair(per):
    """MRR by the number of graphs with evidence for the pair: pooled over
    every KG arm against grounded arms only. A flat grounded line under a
    rising pooled line means the slope is the empty-dossier effect, not pair
    difficulty."""
    d = per.dropna(subset=['n_kg_grounded']).copy()
    d['n_kg_grounded'] = d['n_kg_grounded'].astype(int)
    d = d[d['n_kg_grounded'] > 0]
    pooled = d.groupby('n_kg_grounded')['mrr'].mean()
    grounded = d.groupby('n_kg_grounded')['mrr_grounded'].mean()
    none = per['mrr_none'].mean()
    xs = sorted(d['n_kg_grounded'].unique())
    fig, ax = plt.subplots(figsize=fig_size(1, 1), constrained_layout=True)
    ends = {}
    for name, series, colour, marker, ls in [('Pooled, all KG arms', pooled, NEUTRAL, 'o', '-'),
                                             ('Grounded arms only', grounded, INK, 'D', '-')]:
        ys = [series.get(k, np.nan) for k in xs]
        ax.plot(xs, ys, ls=ls, marker=marker, color=colour, lw=1, ms=3.5, mfc='white',
                mec=colour, mew=0.75, label=name, zorder=3)
        if xs:
            ends[name] = (xs[-1], ys[-1], colour)
    ax.set_xticks(xs)
    ax.set_xticklabels([f'{k} graph' + ('s' if k > 1 else '') for k in xs])
    if xs:
        ax.set_xlim(xs[0] - 0.3, xs[-1] + 0.3)
    ax.set_ylim(0, 1.0)
    ref_line(ax, none, 'no KG')
    clean_ax(ax, xlabel='Graphs with evidence for the pair', ylabel='Mean reciprocal rank on KG arms')
    _end_labels_or_legend(ax, ends)
    return fig


def _odds_rows(df, kept_col, traits):
    kept = df[kept_col].astype(bool)
    rows = [(label, *odds_ratio(mask.astype(bool), kept)) for label, mask in traits.items()]
    rows.sort(key=lambda r: r[1])
    return rows


def fig_bias_forest(df, kept_col, traits, title='Odds of being kept, by trait'):
    """Forest plot of the odds of being kept for every trait in ``traits``
    ({label: boolean Series}), sorted by effect, log scale. A trait whose
    interval excludes one is drawn in ink, the rest in grey. ``title`` is
    accepted for the notebooks and not drawn: the caption carries it."""
    data = _odds_rows(df, kept_col, traits)
    fig, ax = plt.subplots(figsize=(SINGLE_COL_W, 0.18 * len(data) + 0.65), constrained_layout=True)
    for i, (label, orr, lo, hi, p) in enumerate(data):
        sig = p < 0.05
        colour = INK if sig else NEUTRAL
        ax.plot([lo, hi], [i, i], color=colour, lw=0.75, solid_capstyle='butt', zorder=2)
        ax.plot([orr], [i], marker='o', ms=3.5, color=colour, mfc=colour if sig else 'white',
                mec=colour, mew=0.75, zorder=3)
        star = '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'n.s.'
        ax.annotate(star, (hi, i), xytext=(3, 0), textcoords='offset points', va='center',
                    ha='left', fontsize=6, color=colour)
    ax.set_xscale('log')
    # Plain numbers on the log axis: an odds ratio reads as 0.5 or 2, not 2x10^0,
    # and the doubling steps are fixed so a narrow range still gets ticks.
    lo_x, hi_x = ax.get_xlim()
    ax.set_xticks([t for t in (0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32) if lo_x <= t <= hi_x])
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _p: f'{v:g}'))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.minorticks_off()
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels([r[0] for r in data])
    ax.set_ylim(-0.7, len(data) - 0.3)
    _vref(ax, 1, 'OR = 1')
    clean_ax(ax, xlabel='Odds of being kept (95% CI)')
    ax.tick_params(axis='y', length=0)
    return fig


def fig_bias_distribution(df, kept_col, var, title=None):
    """Strip plot of one continuous covariate, kept above dropped, with a bar
    at each median. ``title`` is accepted and not drawn."""
    kept = df[kept_col].astype(bool)
    groups = [('Kept', df[kept], INK), ('Dropped', df[~kept], NEUTRAL)]
    fig, ax = plt.subplots(figsize=(SINGLE_COL_W, 1.3), constrained_layout=True)
    rng = np.random.default_rng(3)
    ticks, tick_labels = [], []
    for i, (label, sub, colour) in enumerate(groups):
        v = sub[var].dropna()
        y0 = 1 - i
        ticks.append(y0)
        tick_labels.append(f'{label} (n = {len(v)})')
        if len(v) == 0:
            continue
        ax.scatter(v, y0 + rng.uniform(-0.18, 0.18, len(v)), s=8, color=colour, alpha=0.5,
                   lw=0, zorder=2)
        med = v.median()
        ax.plot([med, med], [y0 - 0.3, y0 + 0.3], color=colour, lw=1, solid_capstyle='butt',
                zorder=3)
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels)
    ax.set_ylim(-0.6, 1.6)
    # Years are integers; a half-year tick would be a fiction of the locator.
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=8))
    clean_ax(ax, xlabel=human_label(var))
    ax.tick_params(axis='y', length=0)
    return fig


def fig_bias_composition(df, kept_col, cat_vars, title='Composition of kept vs dropped'):
    """Diverging paired bars (dropped to the left, kept to the right) of the
    within-group share of every level, one panel per categorical covariate.
    ``title`` is accepted and not drawn."""
    kept = df[kept_col].astype(bool)
    grp = np.where(kept, 'kept', 'dropped')
    tables = []
    for var in cat_vars:
        ct = pd.crosstab(df[var], grp).reindex(columns=['kept', 'dropped'], fill_value=0)
        ct = ct.loc[ct.sum(axis=1).sort_values().index]
        tables.append((var, ct))
    n_max = max(len(ct) for _, ct in tables)
    fig, axes, ncol = facet_grid(len(tables), row_h=0.16 * n_max + 0.75)
    n_k, n_d = int(kept.sum()), int((~kept).sum())
    # One x scale for every panel, so a bar's length means the same share in each.
    shares = [ct[c] / max(ct[c].sum(), 1) * 100 for _, ct in tables for c in ('kept', 'dropped')]
    mx = max([s.max() for s in shares if len(s)] + [1e-9])
    for ax, (var, ct) in zip(axes, tables):
        k = ct['kept'] / max(ct['kept'].sum(), 1) * 100
        d = ct['dropped'] / max(ct['dropped'].sum(), 1) * 100
        n = len(ct)
        # Panels share a y extent and hang from the top, so the largest
        # level of every covariate sits on the same line across the figure.
        y = np.arange(n) + (n_max - n)
        ax.barh(y, -d.to_numpy(), color=NEUTRAL, height=BAR_W, zorder=2)
        ax.barh(y, k.to_numpy(), color=INK, height=BAR_W, zorder=2)
        ax.axvline(0, color=TEXT_COLOR, lw=0.5, zorder=3)
        # A percentage label is about 12 pt wide; on a 55 mm panel spanning
        # 2.3 times the widest bar it needs a quarter of that bar's length.
        for yi, cat in zip(y, ct.index):
            if d[cat] > 0.25 * mx:
                ax.text(-d[cat] / 2, yi, f'{d[cat]:.0f}%', ha='center', va='center',
                        fontsize=6, color='white')
            if k[cat] > 0.25 * mx:
                ax.text(k[cat] / 2, yi, f'{k[cat]:.0f}%', ha='center', va='center',
                        fontsize=6, color='white')
        ax.set_xlim(-mx * 1.15, mx * 1.15)
        ax.set_ylim(-0.6, n_max - 0.4)
        ax.set_yticks(y)
        ax.set_yticklabels([str(c) for c in ct.index])
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _p: f'{abs(v):.0f}'))
        clean_ax(ax, xlabel='Share of group (%)')
        ax.tick_params(axis='y', length=0)
        facet_title(ax, human_label(var))
    fig.legend(handles=[Patch(color=NEUTRAL, label=f'Dropped (n = {n_d})'),
                        Patch(color=INK, label=f'Kept (n = {n_k})')],
               loc='outside lower center', ncol=2)
    if len(tables) > 1:
        panel_letters(fig, axes)
    return fig


# ───────────────────────── benchmark overview (G01) ─────────────────────────
def fig_benchmark_overview(bench, content_dates, kgs):
    """Three panels a reviewer asks for before any result: when the approvals
    fall relative to each graph's content date, what each graph can evaluate,
    and how many graphs can evaluate each pair.

    ``bench`` is the benchmark frame from ``load_benchmark`` (needs
    ``approval_date``, ``covered_<kg>``, ``leaked_<kg>``, ``grounded_<kg>``,
    ``evaluable_<kg>`` and ``n_kg_grounded``); ``content_dates`` maps a graph to its content date
    (``config/kg_snapshots.yaml``); ``kgs`` are the KG arms in KG_ORDER.
    """
    from biokgsuite.viz.style import KG_NAMES, KG_ORDER
    kgs = [k for k in KG_ORDER if k in kgs]
    fig, axes = plt.subplots(1, 3, figsize=fig_size(2, 1, row_h=2.0), constrained_layout=True,
                             gridspec_kw={'width_ratios': [1.5, 1.0, 0.8]})
    ax_t, ax_s, ax_n = axes

    # a: approvals per quarter, with every graph's content date as a marker.
    # The bars sit in the lower part of the panel; the date labels take the
    # upper part in staggered rows, ordered by date, so a cluster of graphs
    # built in the same year stays legible.
    dates = pd.to_datetime(bench['approval_date'])
    quarters = dates.dt.to_period('Q')
    counts = quarters.value_counts().sort_index()
    x = [p.start_time for p in counts.index]
    ax_t.bar(x, counts.values, width=80, color=NEUTRAL, align='edge')
    ymax = counts.max() * 1.55
    shown = sorted(((k, pd.Timestamp(content_dates[k])) for k in KG_ORDER
                    if content_dates.get(k) is not None), key=lambda t: t[1])
    rows = 3
    for i, (kg, d) in enumerate(shown):
        colour = KG_PALETTE.get(kg, NEUTRAL)
        ax_t.axvline(d, color=colour, lw=0.8, ls='--', zorder=1)
        ax_t.annotate(KG_NAMES.get(kg, kg), xy=(d, ymax), xytext=(2, -1 - 7.5 * (i % rows)),
                      textcoords='offset points', ha='left', va='top', fontsize=6,
                      color=colour, annotation_clip=False)
    lo = min(min(d for _, d in shown), dates.min()) - pd.Timedelta(days=200)
    ax_t.set_xlim(lo, dates.max() + pd.Timedelta(days=150))
    ax_t.set_ylim(0, ymax)
    import matplotlib.dates as mdates
    ax_t.xaxis.set_major_locator(mdates.YearLocator(2))
    ax_t.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    clean_ax(ax_t, xlabel='Approval date (by quarter)', ylabel='Approvals')

    # b: what each graph can evaluate. Evaluable in the graph colour; covered
    # but already carrying the edge (a leak) hatched; covered but with no
    # gene or pathway neighbour to show (not grounded) outlined; not covered
    # FAINT. The legend swatches are neutral because the colour belongs to
    # the graph.
    n = len(bench)
    for i, kg in enumerate(kgs):
        colour = KG_PALETTE.get(kg, NEUTRAL)
        cov, ev = bench[f'covered_{kg}'], bench[f'evaluable_{kg}']
        leaked = bench[f'leaked_{kg}'] if f'leaked_{kg}' in bench.columns else pd.Series(False, index=bench.index)
        grounded = bench[f'grounded_{kg}'] if f'grounded_{kg}' in bench.columns else cov
        n_ev = int(ev.sum())
        n_leak = int((cov & ~ev & leaked).sum())
        n_blank = int((cov & ~ev & ~leaked & ~grounded).sum())
        n_miss = n - n_ev - n_leak - n_blank
        ax_s.barh(i, n_ev, color=colour, height=0.6)
        ax_s.barh(i, n_leak, left=n_ev, color='white', edgecolor=colour, hatch='////', linewidth=0.5, height=0.6)
        ax_s.barh(i, n_blank, left=n_ev + n_leak, color='white', edgecolor=colour, linewidth=0.5, height=0.6)
        ax_s.barh(i, n_miss, left=n_ev + n_leak + n_blank, color=FAINT, height=0.6)
        ax_s.text(n + 2, i, f'{n_ev}', va='center', ha='left', fontsize=6, color=TEXT_COLOR)
    ax_s.set_yticks(range(len(kgs)))
    ax_s.set_yticklabels([KG_NAMES.get(k, k) for k in kgs])
    ax_s.invert_yaxis()
    ax_s.set_xlim(0, n * 1.16)
    clean_ax(ax_s, xlabel='Pairs')
    ax_s.legend(handles=[Patch(color=NEUTRAL, label='evaluable'),
                         Patch(facecolor='white', edgecolor=NEUTRAL, hatch='////', linewidth=0.5,
                               label='covered, edge already present'),
                         Patch(facecolor='white', edgecolor=NEUTRAL, linewidth=0.5,
                               label='covered, no evidence to show'),
                         Patch(color=FAINT, label='not covered')],
                loc='upper center', bbox_to_anchor=(0.5, -0.34), ncol=1, fontsize=6, frameon=False)

    # c: how many graphs have evidence for each pair.
    dist = bench['n_kg_grounded'].value_counts().reindex(range(len(kgs) + 1), fill_value=0)
    ax_n.bar(dist.index, dist.values, width=0.6, color=NEUTRAL)
    ax_n.set_xticks(range(len(kgs) + 1))
    clean_ax(ax_n, xlabel='Graphs with evidence for the pair', ylabel='Pairs')
    panel_letters(fig, axes)
    return fig
