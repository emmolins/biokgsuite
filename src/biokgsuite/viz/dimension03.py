"""Figures for dimension 03, trustworthiness.

One module per dimension, under the package rather than beside the notebook.
Each notebook used to import a bare `figures` module from its own directory via
a sys.path insert, which had two failure modes: the import resolved through the
working directory, so it worked from the notebook's folder and failed from the
repo root; and all five directories named the module `figures`, so Python's
module cache returned whichever was imported first to every notebook that ran
afterwards in the same kernel.
"""
from biokgsuite.viz.primitives import KGS, NICE, facet_grid, kgc, panel_letters, truncate
from biokgsuite.viz.style import clean_ax, facet_title, si_formatter


def fig03_source_diversity(kg_counts, panel_config, top_n=8):
    """Top sources by edge count, one panel per graph.

    Source labels are each graph's own names (Hetionet metaedge codes, BioKG
    predicates) because they are what a reader will look up in that graph's
    documentation; they are only truncated, never translated.
    """
    kgs = [k for k in KGS if k in panel_config and k in kg_counts]
    kgs += [k for k in panel_config if k in kg_counts and k not in kgs]
    if not kgs:
        raise ValueError('no graph in panel_config has counts to plot')

    fig, axes, _ = facet_grid(len(kgs))
    for ax, kg in zip(axes, kgs):
        cfg = panel_config[kg]
        col = cfg['label_col']
        d = kg_counts[kg].sort_values('n_edges', ascending=False).head(top_n).iloc[::-1]
        ax.barh(range(len(d)), d['n_edges'], height=0.6, color=kgc(kg), zorder=2)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels([truncate(s, 28) for s in d[col]])
        ax.set_ylim(-0.6, len(d) - 0.4)
        ax.set_xlim(0, d['n_edges'].max() * 1.05)
        ax.xaxis.set_major_formatter(si_formatter())
        ax.locator_params(axis='x', nbins=4)
        clean_ax(ax, xlabel='Edge count')
        facet_title(ax, NICE.get(kg, cfg.get('title', kg)))

    if len(kgs) > 1:
        panel_letters(fig, axes)
    return fig
