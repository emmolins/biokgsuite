"""Figures for dimension 03, trustworthiness.

One module per dimension, under the package rather than beside the notebook.
Each notebook used to import a bare `figures` module from its own directory via
a sys.path insert, which had two failure modes: the import resolved through the
working directory, so it worked from the notebook's folder and failed from the
repo root; and all five directories named the module `figures`, so Python's
module cache returned whichever was imported first to every notebook that ran
afterwards in the same kernel.
"""
import math

import numpy as np
import matplotlib.pyplot as plt

from biokgsuite.viz.style import TEXT_COLOR
from biokgsuite.viz.primitives import kgc


def fig03_source_diversity(kg_counts, panel_config, top_n=8):
    """Top sources by edge count, one panel per graph."""
    kgs = [k for k in panel_config if k in kg_counts]
    if not kgs:
        raise ValueError('no graph in panel_config has counts to plot')
    ncol = min(3, len(kgs))
    nrow = math.ceil(len(kgs) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.2 * nrow + 0.5))
    axes = np.atleast_1d(np.array(axes)).reshape(-1)

    for ax, kg in zip(axes, kgs):
        cfg = panel_config[kg]
        col = cfg['label_col']
        d = kg_counts[kg].sort_values('n_edges', ascending=False).head(top_n).iloc[::-1]
        lab = [str(s)[:24] + ('…' if len(str(s)) > 24 else '') for s in d[col]]
        ax.barh(range(len(d)), d['n_edges'], color=kgc(kg),
                edgecolor='white', lw=0.5, zorder=2)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels(lab, fontsize=7.5)
        ax.set_xlabel('edge count', fontsize=8.5)
        ax.tick_params(axis='x', labelsize=7.5)
        ax.set_title(cfg['title'], fontsize=10.5, fontweight='bold',
                     loc='left', color=TEXT_COLOR, pad=6)
        ax.set_xlim(0, d['n_edges'].max() * 1.18)

    for ax in axes[len(kgs):]:
        ax.set_visible(False)

    fig.suptitle('Source-database diversity across knowledge graphs', fontsize=13,
                 fontweight='bold', x=0.015, ha='left', color=TEXT_COLOR, y=0.99)
    fig.text(0.015, 0.945, f'top {top_n} sources by edge count per knowledge graph',
             fontsize=9.5, color='#777', ha='left')
    fig.subplots_adjust(left=0.13, right=0.97, top=0.88, bottom=0.09,
                        wspace=0.55, hspace=0.55)
    return fig
