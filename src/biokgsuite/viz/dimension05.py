"""Figures for dimension 05, stability.

One module per dimension, under the package rather than beside the notebook, so
there is a single import root and no `sys.path` manipulation. Each function
returns its figure rather than saving it; saving belongs to
`biokgsuite.io.Outputs`, so every notebook writes through one path, uses one
naming convention, and reports what it produced.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from biokgsuite.viz.primitives import NICE, _present, facet_grid, kgc, panel_letters, ref_line
from biokgsuite.viz.style import NEUTRAL, TEXT_COLOR, clean_ax, facet_title, fig_size

# Interpretive reference lines, not scoring thresholds.
HIGH_STABILITY = 0.95
FRAGILE = 0.85


def fig05_stability(raw_r, dropout_rates, raw_r_null=None):
    """Spearman r vs edge dropout, one panel per KG.

    Two encoding channels, held constant across every panel:
      colour     the pair set   (KG colour for therapeutic pairs, grey for null)
      linestyle  the strategy   (solid for random dropout, dashed for periphery)

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

    fig, axes, ncol = facet_grid(len(kgs))

    def series(ax, ys, color, ls, marker):
        ax.plot(pct, ys, ls=ls, marker=marker, color=color, lw=1, ms=3.5,
                markeredgecolor='white', markeredgewidth=0.5, zorder=3)

    for idx, (ax, kg) in enumerate(zip(axes, kgs)):
        c = kgc(kg)
        series(ax, [raw_r[kg]['random'].get(r, np.nan) for r in rates], c, '-', 'o')
        series(ax, [raw_r[kg]['periphery'].get(r, np.nan) for r in rates], c, '--', 's')
        if raw_r_null and kg in raw_r_null:
            series(ax, [raw_r_null[kg]['random'].get(r, np.nan) for r in rates],
                   NEUTRAL, '-', 'o')
            series(ax, [raw_r_null[kg]['periphery'].get(r, np.nan) for r in rates],
                   NEUTRAL, '--', 's')
        ref_line(ax, HIGH_STABILITY, f'r = {HIGH_STABILITY}')
        ref_line(ax, FRAGILE, f'r = {FRAGILE}')
        ax.set_xticks(pct)
        ax.set_xticklabels([f'{p:g}' for p in pct])
        ax.set_ylim(min(lo - 0.02, FRAGILE - 0.02), 1.005)
        clean_ax(ax, xlabel='Edge dropout (%)' if idx >= len(kgs) - ncol else '',
                 ylabel='Spearman r' if idx % ncol == 0 else '')
        if idx % ncol:
            ax.tick_params(labelleft=False)
        facet_title(ax, NICE.get(kg, kg))

    # One legend for the whole figure, below the panels, so it never sits over
    # data and is not repeated. The therapeutic swatches are drawn in neutral
    # ink because their colour is the graph's, which the facet title names.
    handles = [
        Line2D([0], [0], color=TEXT_COLOR, lw=1, ls='-', marker='o', ms=3.5,
               label='Random dropout'),
        Line2D([0], [0], color=TEXT_COLOR, lw=1, ls='--', marker='s', ms=3.5,
               label='Periphery dropout'),
    ]
    if raw_r_null:
        handles.append(Line2D([0], [0], color=NEUTRAL, lw=1, ls='-', marker='o', ms=3.5,
                              label='Degree-matched null'))
    fig.legend(handles=handles, loc='outside lower center', ncol=len(handles))
    if len(kgs) > 1:
        panel_letters(fig, axes)
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

    fig, axes = plt.subplots(1, len(strategies), figsize=fig_size(2, 1),
                             sharey=True, squeeze=False, constrained_layout=True)
    axes = axes.reshape(-1)
    # Each scorer owns one unit slot; the graphs share 0.8 of it, and within a
    # graph the therapeutic bar and its null bar sit side by side.
    slot = 0.8 / max(1, len(kgs))
    width = slot * 0.45

    for ax, strategy in zip(axes, strategies):
        sub = df[df['Strategy'] == strategy]
        x = np.arange(len(scorer_names))
        for k, kg in enumerate(kgs):
            centre = x + (k - (len(kgs) - 1) / 2) * slot
            for j, (pair_set, colour) in enumerate((('gold', kgc(kg)), ('null', NEUTRAL))):
                vals = [sub[(sub['KG'] == kg) & (sub['Scorer'] == sc)
                            & (sub['Pair set'] == pair_set)]['Spearman r'].mean()
                        for sc in scorer_names]
                ax.bar(centre + (j - 0.5) * width, vals, width, color=colour, zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels([sc + ('\n(scored)' if sc == scored_by else '')
                            for sc in scorer_names])
        ax.set_ylim(min(0.8, float(df['Spearman r'].min()) - 0.02), 1.005)
        clean_ax(ax, ylabel='Spearman r' if ax is axes[0] else '')
        facet_title(ax, strategy)

    handles = ([Patch(color=kgc(kg), label=NICE.get(kg, kg)) for kg in kgs]
               + [Patch(color=NEUTRAL, label='Degree-matched null')])
    fig.legend(handles=handles, loc='outside lower center', ncol=len(handles),
               handlelength=1.0)
    panel_letters(fig, axes)
    return fig
