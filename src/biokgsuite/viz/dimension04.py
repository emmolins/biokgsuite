"""Figures for dimension 04, topology.

One module per dimension, under the package rather than beside the notebook, so
there is a single import root and no `sys.path` manipulation. All five dimension
folders previously defined a module named `figures`, and Python's module cache
returns whichever was imported first to every notebook that runs afterwards in
the same kernel.

Each function RETURNS its figure rather than saving it. Saving belongs to
`biokgsuite.io.Outputs`, so every notebook writes through one path, uses one
naming convention, and reports what it produced. Functions that saved
themselves bypassed that accounting entirely: notebook 04's four figures were
silently absent from the manifest because nothing tracked them.
"""
import numpy as np
import matplotlib.pyplot as plt

from biokgsuite.viz.primitives import (NICE, _df, _present, _xlabels, facet_grid, kgc,
                                       panel_letters, ref_line)
from biokgsuite.viz.style import (MUTED, NEUTRAL, TEXT_COLOR, clean_ax, facet_title,
                                  fig_size)

BAR_W = 0.6


def _values_above(ax, xs, vals, texts, log=False):
    """6 pt value above each bar; used only where the number is the message."""
    for x, v, t in zip(xs, vals, texts):
        y = v * 1.12 if log else v + 0.012 * ax.get_ylim()[1]
        ax.text(x, y, t, ha='center', va='bottom', fontsize=6, color=TEXT_COLOR)


# ───────────────────────── 04 · clustering coefficient ──────────────────────────
def fig04_clustering(clust):
    """Clustering coefficient relative to an Erdos-Renyi graph of the same size.

    The multiplier is the finding (57x, 380x), so it is printed above each bar;
    the log axis alone would hide a five-fold difference between graphs.
    """
    cl = _df(clust)
    cl = cl.reindex(_present(cl.index))
    order = list(cl.index)
    vals = cl['C / C_ER'].values
    xs = np.arange(len(order))

    fig, ax = plt.subplots(figsize=fig_size(1, 1), constrained_layout=True)
    ax.bar(xs, vals, width=BAR_W, color=[kgc(k) for k in order], zorder=2)
    ax.set_yscale('log')
    # Head room for the value labels, then the 1x line labelled past the axis.
    ax.set_ylim(min(0.5, np.nanmin(vals) / 2), np.nanmax(vals) * 3)
    ax.set_xlim(-0.6, len(order) - 0.4)
    _values_above(ax, xs, vals, [f'{v:.0f}×' for v in vals], log=True)
    ref_line(ax, 1, '1× (random)')
    _xlabels(ax, order)
    clean_ax(ax, ylabel='Clustering vs random, C / C$_{ER}$')
    return fig


# ───────────────────────── 04 · graph cohesion ──────────────────────────────────
def fig04_cohesion(conn):
    """Share of nodes in the largest component, and the number of components."""
    conn = _df(conn)
    conn = conn.reindex(_present(conn.index))
    order = list(conn.index)
    cols = [kgc(k) for k in order]
    xs = np.arange(len(order))

    fig, axes = plt.subplots(1, 2, figsize=fig_size(2, 1), constrained_layout=True)

    ax = axes[0]
    lcc = conn['LCC (%)'].values
    ax.bar(xs, lcc, width=BAR_W, color=cols, zorder=2)
    ax.set_ylim(0, 112)
    ax.set_xlim(-0.6, len(order) - 0.4)
    _values_above(ax, xs, lcc, [f'{v:.1f}' for v in lcc])
    ref_line(ax, 95, '95%')
    _xlabels(ax, order)
    clean_ax(ax, ylabel='Nodes in largest component (%)')
    facet_title(ax, 'Largest connected component')

    ax = axes[1]
    comp = conn['Components'].values
    ax.bar(xs, comp, width=BAR_W, color=cols, zorder=2)
    ax.set_yscale('log')
    ax.set_ylim(1, np.nanmax(comp) * 4)
    ax.set_xlim(-0.6, len(order) - 0.4)
    _values_above(ax, xs, comp, [f'{int(v):,}' for v in comp], log=True)
    _xlabels(ax, order)
    clean_ax(ax, ylabel='Connected components')
    facet_title(ax, 'Connected components')

    panel_letters(fig, axes)
    return fig


# ───────────────────────── 04c · known-pair recovery ────────────────────────────
def fig04_recovery(recovery_results):
    """Share of known drug-disease pairs reachable within k hops, per graph.

    Graphs are named at the line ends when the ends are far enough apart to
    read; otherwise (most graphs saturate near 100% at k = 5) a legend.
    """
    ks = [1, 2, 3, 4, 5]
    kgs = _present(recovery_results.keys())
    fig, ax = plt.subplots(figsize=fig_size(1, 1), constrained_layout=True)
    ends = {}
    for kg in kgs:
        r = recovery_results.get(kg, {}).get('recovery', {})
        ys = [r.get(k, np.nan) * 100 for k in ks]
        ax.plot(ks, ys, marker='o', ms=3.5, lw=1, color=kgc(kg), mfc='white',
                mec=kgc(kg), mew=0.75, label=NICE.get(kg, kg), zorder=3)
        ends[kg] = ys[-1]
    ax.set_xticks(ks)
    ax.set_ylim(0, 105)
    ax.set_xlim(0.8, 5.2)
    clean_ax(ax, xlabel='Hop distance (k)', ylabel='Known pairs recovered (%)')

    # Direct labels need one text height between neighbouring line ends;
    # 6 pt text is about 6% of a 50 mm panel, so ask for 7% of the y range.
    finite = sorted(v for v in ends.values() if v == v)
    gaps = np.diff(finite) if len(finite) > 1 else np.array([np.inf])
    if kgs and gaps.min() >= 0.07 * 105:
        ax.set_xlim(0.8, 5.9)
        for kg, y in ends.items():
            if y == y:
                ax.text(5.12, y, NICE.get(kg, kg), fontsize=6, color=TEXT_COLOR,
                        ha='left', va='center')
    elif len(kgs) > 1:
        ax.legend(loc='lower right', ncol=1 if len(kgs) <= 3 else 2)
    return fig


# ───────────────────────── 04 · differential resilience (dropout curves) ────────
def fig04_resilience(diversity_results):
    """Recovery of gold-standard versus random pairs as edges are dropped.

    The band between the two curves is the resilience the graph owes to
    therapeutic structure rather than to density; DR summarises its area.
    """
    from matplotlib.lines import Line2D
    kgs = _present(diversity_results.keys())
    fig, axes, ncol = facet_grid(len(kgs))
    for idx, (ax, kg) in enumerate(zip(axes, kgs)):
        r = diversity_results[kg]
        x = [v * 100 for v in r['drop_rates']]
        gc, rc = list(r['gold_curve']), list(r['rand_curve'])
        c = kgc(kg)
        ax.plot(x, gc, '-o', color=c, lw=1, ms=3.5, zorder=3)
        ax.plot(x, rc, '--s', color=NEUTRAL, lw=1, ms=3.5, zorder=3)
        ax.fill_between(x, rc, gc, where=[g >= rr for g, rr in zip(gc, rc)],
                        alpha=0.12, color=c, lw=0, zorder=2)
        ax.set_ylim(0, 1.05)
        ax.set_xlim(-2, max(x) + 2)
        ax.set_xticks(x)
        # The score sits bottom right, lifted just above whichever curve
        # reaches into that corner so it never lands on a marker.
        right = [i for i, xv in enumerate(x) if xv >= max(x) * 0.7]
        floor = max([rc[i] for i in right] + [gc[i] for i in right if gc[i] < 0.3] + [0])
        ax.text(0.97, max(0.04, floor / 1.05 + 0.04), f"DR {r['dr_score']:.2f}",
                transform=ax.transAxes, ha='right', va='bottom', fontsize=6, color=MUTED)
        clean_ax(ax, xlabel='Edge dropout (%)' if idx >= len(kgs) - ncol else '',
                 ylabel='Recovery fraction' if idx % ncol == 0 else '')
        if idx % ncol:
            ax.tick_params(labelleft=False)
        facet_title(ax, NICE.get(kg, kg))

    # The gold series is a different colour in every facet, so its legend
    # swatch is neutral ink; the facet title already says which graph it is.
    handles = [Line2D([0], [0], color=TEXT_COLOR, lw=1, marker='o', ms=3.5,
                      label='Gold-standard pairs'),
               Line2D([0], [0], color=NEUTRAL, lw=1, ls='--', marker='s', ms=3.5,
                      label='Random pairs')]
    fig.legend(handles=handles, loc='outside lower center', ncol=2)
    if len(kgs) > 1:
        panel_letters(fig, axes)
    return fig
