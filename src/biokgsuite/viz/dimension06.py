"""Figures for dimension 06, task performance.

One module per dimension, under the package rather than beside the notebook, so
there is a single import root and no `sys.path` manipulation. All the dimension
folders previously defined a module named `figures`, and Python's module cache
returns whichever was imported first to every notebook that runs afterwards in
the same kernel.

Each function RETURNS its figure rather than saving it. Saving belongs to
`biokgsuite.io.Outputs`, so every notebook writes through one path and reports
what it produced.

Three of these four figures used to live inline in the notebook and save
themselves, which is why none of them appeared in the results manifest and why
none could be redrawn without re-running the analysis. The precision-recall
panel additionally RE-SCORED every pair with every heuristic inside the figure
cell, duplicating the most expensive work in the notebook and leaving the drawn
curves free to drift from the AUPRC they were labelled with. It now takes the
curves computed from the same stored scores the metrics came from.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from biokgsuite.viz.style import (clean_ax, panel_label, HEATMAP_CMAP,
                                  HEURISTIC_COLORS, TEXT_COLOR)
from biokgsuite.viz.primitives import kgc, _present, NICE

STRATEGY_LABEL = {'random': 'Random', 'type-constrained': 'Type-constrained',
                  'shared-target': 'Shared-target'}
STRATEGY_ORDER = ['random', 'type-constrained', 'shared-target']


def _left_spine(ax):
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_color('#333333')
    ax.spines['left'].set_linewidth(0.5)


# ───────────────────────── 06a · link prediction, heatmap ───────────────────────
def fig06_link_heatmap(link_records, primary_heuristic):
    """AUROC and AUPRC per KG per negative-sampling strategy, primary heuristic."""
    df = pd.DataFrame(link_records)
    df = df[df['heuristic'] == primary_heuristic]
    strat = [s for s in STRATEGY_ORDER if s in set(df['strategy'])] or list(df['strategy'].unique())
    kgs = _present(df['kg'].unique())

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6))
    for ax, metric in zip(axes, ['auroc', 'auprc']):
        piv = (df.pivot_table(index='kg', columns='strategy', values=metric, aggfunc='mean')
                 .reindex(index=kgs, columns=strat))
        im = ax.pcolormesh(piv.values.astype(float), cmap=HEATMAP_CMAP,
                           vmin=0.45, vmax=1.0, edgecolors='white', linewidth=2)
        ax.set_xticks([x + 0.5 for x in range(len(strat))])
        ax.set_xticklabels([STRATEGY_LABEL.get(s, s) for s in strat], rotation=20, ha='right')
        ax.set_yticks([y + 0.5 for y in range(len(kgs))])
        ax.set_yticklabels([NICE.get(k, k) for k in kgs])
        ax.tick_params(length=0); ax.invert_yaxis(); ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        for i in range(len(kgs)):
            for j in range(len(strat)):
                v = piv.iloc[i, j]
                if pd.notna(v):
                    ax.text(j + 0.5, i + 0.5, f'{v:.3f}', ha='center', va='center',
                            fontsize=9, fontweight='bold',
                            color='white' if v > 0.7 else TEXT_COLOR)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(f'{metric.upper()}  ({primary_heuristic})', fontsize=10.5,
                     fontweight='bold', loc='left', color=TEXT_COLOR, pad=8)
    fig.suptitle('Link-prediction performance by negative-sampling strategy',
                 fontsize=12, fontweight='bold', x=0.012, ha='left',
                 color=TEXT_COLOR, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


# ───────────────────────── 06b · link prediction, PR curves ─────────────────────
def fig06_link_pr(pr_curves, strategy_label='type-constrained'):
    """Precision-recall per KG, one line per heuristic, against prevalence.

    `pr_curves` maps kg -> {'baseline': float,
                            'curves': {heuristic: (recall, precision, auprc)}}
    built from the SAME scores the AUPRC was computed from, so a curve and its
    label cannot disagree.
    """
    kgs = _present(pr_curves.keys())
    fig, axes = plt.subplots(1, max(1, len(kgs)), figsize=(4.3 * max(1, len(kgs)), 3.8),
                             sharey=True, squeeze=False)
    axes = axes.reshape(-1)

    for ax, name in zip(axes, kgs):
        entry = pr_curves[name]
        for heuristic, (rec, prec, _auprc) in entry['curves'].items():
            ax.plot(rec, prec, color=HEURISTIC_COLORS.get(heuristic, '#9e9288'),
                    lw=1.8, label=heuristic)
        ax.axhline(entry['baseline'], color='#868074', ls='--', lw=1.2, alpha=0.85,
                   label='Prevalence')
        # AUPRC and prevalence are PER PANEL and belong in the panel. The shared
        # legend used to carry the first KG's numbers while sitting under every
        # panel, so Hetionet's AUPRC read as if it described BioKG too.
        lines = [f'{h}  {a:.3f}' for h, (_, _, a) in entry['curves'].items()]
        lines.append(f"prevalence  {entry['baseline']:.3f}")
        # Lower-left: a precision-recall curve is high at low recall and falls at
        # high recall, so this corner is the one reliably free of data.
        ax.text(0.03, 0.03, 'AUPRC\n' + '\n'.join(lines), transform=ax.transAxes,
                ha='left', va='bottom', fontsize=7.5, color=TEXT_COLOR,
                linespacing=1.45,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          edgecolor='#dddddd', linewidth=0.6, alpha=0.9))
        ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.06)
        _left_spine(ax)
        clean_ax(ax, title=NICE.get(name, name.upper()), xlabel='Recall',
                 ylabel='Precision' if ax is axes[0] else '')

    for ax in axes[len(kgs):]:
        ax.set_visible(False)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.06),
               ncol=4, fontsize=9, frameon=False, handlelength=1.6,
               handletextpad=0.5, columnspacing=1.4)
    fig.suptitle(f'Link prediction: precision-recall ({strategy_label} negatives)',
                 fontsize=11, fontweight='bold', y=1.04, color=TEXT_COLOR)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.20)
    return fig


# ───────────────────────── 06c · neighbourhood retrieval ────────────────────────
def fig06_retrieval(target_results, k_range=(1, 5, 10, 20, 50, 100, 200, 500)):
    """Recall@K curves, and Recall@10/50/100 as grouped bars."""
    kgs = _present([n for n in target_results if 'ranks' in target_results[n]])
    k_range = list(k_range)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    for name in kgs:
        ranks = np.asarray(target_results[name]['ranks'])
        ax.plot(k_range, [float(np.mean(ranks <= k)) for k in k_range], 'o-',
                color=kgc(name), lw=2.2, markersize=5.5, label=NICE.get(name, name.upper()))
    clean_ax(ax, title='Recall@K', xlabel='K', ylabel='Recall')
    ax.set_xlim(0, max(k_range) * 1.05); ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8, frameon=True, facecolor='white', edgecolor='none',
              framealpha=0.9, loc='upper left', bbox_to_anchor=(0.02, 0.98))
    _left_spine(ax); panel_label(ax, 'a')

    ax = axes[1]
    ks = ['recall@10', 'recall@50', 'recall@100']
    x, w = np.arange(len(ks)), 0.15
    for i, name in enumerate(kgs):
        vals = [target_results[name][k] for k in ks]
        bars = ax.bar(x + i * w, vals, w, label=NICE.get(name, name.upper()),
                      color=kgc(name), edgecolor='none', linewidth=0)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                    f'{v:.2f}', ha='center', fontsize=8, color=TEXT_COLOR,
                    rotation=90 if len(kgs) > 4 else 0)
    ax.set_xticks(x + w * (len(kgs) - 1) / 2)
    ax.set_xticklabels(['Recall@10', 'Recall@50', 'Recall@100'], fontsize=10)
    ax.set_ylim(0, 1.15)
    clean_ax(ax, title='Recall summary', ylabel='Recall')
    _left_spine(ax); panel_label(ax, 'b')

    fig.suptitle('Neighbourhood retrieval', fontsize=12, fontweight='bold',
                 y=1.02, color=TEXT_COLOR)
    fig.tight_layout()
    return fig


# ───────────────────────── 06d · multi-hop reasoning ────────────────────────────
def fig06_multihop(multihop_results):
    """Hits@K for the intermediate target, and the rank CDF behind it."""
    kgs = _present([n for n, r in multihop_results.items()
                    if r.get('n_evaluated', 0) > 0])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    ks = ['hits@10', 'hits@50', 'hits@100']
    x, w = np.arange(len(ks)), 0.15
    for i, name in enumerate(kgs):
        vals = [multihop_results[name][k] for k in ks]
        bars = ax.bar(x + i * w, vals, w, label=NICE.get(name, name.upper()),
                      color=kgc(name), edgecolor='none', linewidth=0)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                    f'{v:.2f}', ha='center', fontsize=8.5, color=TEXT_COLOR, rotation=90)
    ax.set_xticks(x + w * (len(kgs) - 1) / 2)
    ax.set_xticklabels(['Hits@10', 'Hits@50', 'Hits@100'], fontsize=10)
    ax.set_ylim(0, 1.15)
    clean_ax(ax, title='Intermediate target recovery', ylabel='Fraction recovered')
    ax.legend(fontsize=9, frameon=False, loc='upper left')
    _left_spine(ax); panel_label(ax, 'a')

    ax = axes[1]
    for name in kgs:
        r = np.sort(np.asarray(multihop_results[name]['ranks']))
        ax.plot(r, np.arange(1, len(r) + 1) / len(r), '-', color=kgc(name), lw=2.2,
                label=f"{NICE.get(name, name.upper())} (MRR={multihop_results[name]['mrr']:.3f})")
    ax.set_xscale('log')
    ax.set_xlabel('Rank of known target (log scale)', fontsize=10)
    ax.set_ylabel('Cumulative fraction', fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.axvline(100, color='#868074', ls='--', lw=1.0, alpha=0.6, label='K = 100')
    clean_ax(ax, title='Target rank CDF')
    ax.legend(fontsize=8, frameon=True, facecolor='white', edgecolor='none',
              framealpha=0.9, loc='upper left', bbox_to_anchor=(0.02, 0.98))
    _left_spine(ax); panel_label(ax, 'b')

    fig.suptitle('Multi-hop reasoning: intermediate target identification',
                 fontsize=12, fontweight='bold', y=1.02, color=TEXT_COLOR)
    fig.tight_layout()
    return fig
