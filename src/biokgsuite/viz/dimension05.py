"""Figures for dimension 05, stability.

One module per dimension, under the package rather than beside the notebook, so
there is a single import root and no `sys.path` manipulation. All the dimension
folders previously defined a module named `figures`, and Python's module cache
returns whichever was imported first to every notebook that runs afterwards in
the same kernel.

Each function RETURNS its figure rather than saving it. Saving belongs to
`biokgsuite.io.Outputs`, so every notebook writes through one path, uses one
naming convention, and reports what it produced.
"""
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from biokgsuite.viz.style import TEXT_COLOR
from biokgsuite.viz.primitives import kgc, _present, NICE

# Interpretive reference lines, not scoring thresholds.
HIGH_STABILITY = 0.95
FRAGILE = 0.85

# The null is drawn in grey in every panel: colour carries the pair set, not the KG.
NULL_COLOR = '#A8A29A'


def fig05_stability(raw_r, dropout_rates, raw_r_null=None):
    """Spearman r vs edge dropout, one panel per KG.

    Two encoding channels, held constant across every panel:
      colour     the pair set   -- KG colour for therapeutic pairs, grey for null
      linestyle  the strategy   -- solid for random dropout, dashed for periphery

    `raw_r_null` draws the degree-matched null. Adamic-Adar is degree-driven, so
    the therapeutic curve alone cannot be read as evidence of robust therapeutic
    structure; the gap between a curve and its null is what the figure is for,
    and the gap changes sign between the two strategies.
    """
    kgs = _present(raw_r.keys())
    rates = list(dropout_rates)
    pct = [r * 100 for r in rates]

    sources = [raw_r] + ([raw_r_null] if raw_r_null else [])
    values = [src[kg][s][r] for src in sources for kg in kgs
              for s in ('random', 'periphery')
              for r in rates if r in src.get(kg, {}).get(s, {})]
    finite = [v for v in values if v == v]
    lo = min([*finite, FRAGILE]) if finite else FRAGILE

    ncol = min(3, max(1, len(kgs)))
    nrow = math.ceil(len(kgs) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.7 * nrow),
                             sharey=True, squeeze=False)
    axes = axes.reshape(-1)

    def series(ax, ys, color, ls, marker):
        ax.plot(pct, ys, ls=ls, marker=marker, color=color, lw=1.9, ms=5,
                markeredgecolor='white', markeredgewidth=0.6, zorder=3)

    for idx, (ax, kg) in enumerate(zip(axes, kgs)):
        c = kgc(kg)
        series(ax, [raw_r[kg]['random'].get(r, np.nan) for r in rates], c, '-', 'o')
        series(ax, [raw_r[kg]['periphery'].get(r, np.nan) for r in rates], c, '--', 's')
        if raw_r_null and kg in raw_r_null:
            series(ax, [raw_r_null[kg]['random'].get(r, np.nan) for r in rates],
                   NULL_COLOR, '-', 'o')
            series(ax, [raw_r_null[kg]['periphery'].get(r, np.nan) for r in rates],
                   NULL_COLOR, '--', 's')
        ax.axhline(HIGH_STABILITY, color='#655F55', ls=':', lw=0.9, alpha=0.8, zorder=1)
        ax.axhline(FRAGILE, color='#E69F00', ls=':', lw=0.9, alpha=0.8, zorder=1)
        ax.set_xticks(pct)
        ax.set_xticklabels([f'{p:g}%' for p in pct])
        ax.set_ylim(min(lo - 0.02, FRAGILE - 0.02), 1.005)
        if idx % ncol == 0:
            ax.set_ylabel('Spearman r')
        if idx >= len(kgs) - ncol:
            ax.set_xlabel('edge dropout (%)')
        ax.set_title(NICE.get(kg, kg), fontsize=10.5, fontweight='bold',
                     loc='left', color=TEXT_COLOR, pad=6)

    for ax in axes[len(kgs):]:
        ax.set_visible(False)

    # One legend for the whole figure, below the panels, so it never sits over
    # data and does not have to be repeated. Colour and linestyle are named
    # separately because they are independent channels.
    # The legend carries only the linestyle channel. The colour channel is stated
    # in the subtitle instead: a legend swatch for "therapeutic pairs" would have
    # to be drawn in some single colour, and every such colour either collides
    # with the null's grey or misrepresents a per-KG palette.
    swatch = '#555555'
    handles = [
        Line2D([0], [0], color=swatch, lw=1.9, ls='-', marker='o', ms=5,
               label='random dropout'),
        Line2D([0], [0], color=swatch, lw=1.9, ls='--', marker='s', ms=5,
               label='periphery dropout'),
        Line2D([0], [0], color='#655F55', lw=0.9, ls=':', label=f'r = {HIGH_STABILITY}'),
        Line2D([0], [0], color='#E69F00', lw=0.9, ls=':', label=f'r = {FRAGILE}'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.005))

    fig.suptitle('Ranking stability under edge dropout', fontsize=12.5,
                 fontweight='bold', x=0.008, ha='left', color=TEXT_COLOR, y=0.985)
    if raw_r_null:
        fig.text(0.008, 0.905,
                 'coloured: therapeutic pairs      grey: degree-matched null '
                 '(same endpoint degrees, no therapeutic relation)',
                 ha='left', fontsize=9, color=TEXT_COLOR, alpha=0.8)
    fig.subplots_adjust(left=0.085, right=0.98, top=0.84, bottom=0.22,
                        wspace=0.12, hspace=0.45)
    return fig


def fig05_scorers(stab_df, primary_rate, scorer_names, scored_by):
    """Do the three neighbourhood scorers agree, or is this an Adamic-Adar fact?

    One panel per strategy; therapeutic pairs and their degree-matched null side
    by side at the scored rate. If the bars move together the conclusion is a
    property of the graph; if Adamic-Adar stands apart, it is a property of its
    inverse-log-degree weighting and the manuscript has to say so.
    """
    df = stab_df[stab_df['Dropout'] == primary_rate]
    kgs = _present(df['KG'].unique())
    strategies = list(dict.fromkeys(df['Strategy']))

    fig, axes = plt.subplots(1, len(strategies),
                             figsize=(5.0 * len(strategies), 3.9),
                             sharey=True, squeeze=False)
    axes = axes.reshape(-1)
    width = 0.36

    for ax, strategy in zip(axes, strategies):
        sub = df[df['Strategy'] == strategy]
        x = np.arange(len(scorer_names))
        for k, kg in enumerate(kgs):
            for j, (pair_set, alpha, hatch) in enumerate(
                    (('gold', 1.0, None), ('null', 0.45, '///'))):
                vals = [sub[(sub['KG'] == kg) & (sub['Scorer'] == sc)
                            & (sub['Pair set'] == pair_set)]['Spearman r'].mean()
                        for sc in scorer_names]
                ax.bar(x + (k - (len(kgs) - 1) / 2) * width + (j - 0.5) * width / 2,
                       vals, width / 2, color=kgc(kg), alpha=alpha, hatch=hatch,
                       edgecolor='white', linewidth=0.5, zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels([sc + ('\n(scored)' if sc == scored_by else '')
                            for sc in scorer_names], fontsize=8.5)
        ax.set_title(strategy, fontsize=10.5, fontweight='bold', loc='left',
                     color=TEXT_COLOR, pad=6)
        ax.set_ylim(min(0.8, float(df['Spearman r'].min()) - 0.02), 1.005)
    axes[0].set_ylabel('Spearman r')

    handles = ([Line2D([0], [0], color=kgc(kg), lw=6, label=NICE.get(kg, kg))
                for kg in kgs]
               + [Line2D([0], [0], color='#555555', lw=6, label='therapeutic pairs'),
                  Line2D([0], [0], color='#555555', lw=6, alpha=0.45,
                         label='degree-matched null (hatched)')])
    fig.legend(handles=handles, loc='lower center', ncol=len(handles), fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(f'Scorer agreement at {primary_rate:.0%} dropout', fontsize=12.5,
                 fontweight='bold', x=0.008, ha='left', color=TEXT_COLOR, y=0.985)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.85, bottom=0.26, wspace=0.10)
    return fig
