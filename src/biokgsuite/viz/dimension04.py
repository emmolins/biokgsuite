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

from biokgsuite.viz.style import TEXT_COLOR, ALERT_RED
from biokgsuite.viz.primitives import kgc, _present, _xlabels, _title, _df, NICE

# ───────────────────────── 04 · clustering coefficient ──────────────────────────
def fig04_clustering(clust):
    cl = _df(clust); cl = cl.reindex(_present(cl.index)); order = list(cl.index)
    cols = [kgc(k) for k in order]
    fig, ax = plt.subplots(figsize=(7.4, 4.7))
    ax.bar(range(len(order)), cl['C / C_ER'], color=cols, edgecolor='white', linewidth=0.5, zorder=2); ax.set_yscale('log')
    for i, val in enumerate(cl['C / C_ER']):
        ax.text(i, val * 1.18, f'{val:.0f}×', ha='center', fontsize=8, fontweight='bold', color=TEXT_COLOR)
    ax.axhline(1, color='#bbb', ls=':', lw=1)
    _xlabels(ax, order); ax.set_ylabel('clustering vs random  (C / C$_{ER}$, log)')
    _title(ax, 'Clustering coefficient relative to a random graph', 'values > 1 indicate small-world structure', pad=12)
    fig.tight_layout(); return fig


# ───────────────────────── 04c · known-pair recovery ────────────────────────────


# ───────────────────────── 04 · graph cohesion ──────────────────────────────────
def fig04_cohesion(conn):
    conn = _df(conn); conn = conn.reindex(_present(conn.index)); order = list(conn.index)
    cols = [kgc(k) for k in order]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6))
    ax = axes[0]
    ax.bar(range(len(order)), conn['LCC (%)'], color=cols, edgecolor='white', linewidth=0.5, zorder=2)
    ax.axhline(95, color=ALERT_RED, ls='--', lw=1.2, label='95% threshold'); ax.set_ylim(0, 112)
    for i, val in enumerate(conn['LCC (%)']):
        ax.text(i, val + 1.5, f'{val:.1f}%', ha='center', fontsize=8, fontweight='bold', color=TEXT_COLOR)
    _xlabels(ax, order); ax.set_ylabel('% of nodes in LCC'); ax.legend(loc='lower left')
    _title(ax, 'Largest connected component')
    ax = axes[1]
    comp = conn['Components'].values
    ax.bar(range(len(order)), comp, color=cols, edgecolor='white', linewidth=0.5, zorder=2); ax.set_yscale('log')
    for i, val in enumerate(comp):
        ax.text(i, val * 1.25, f'{int(val):,}', ha='center', fontsize=8, fontweight='bold', color=TEXT_COLOR)
    _xlabels(ax, order); ax.set_ylabel('component count (log scale)')
    _title(ax, 'Connected components')
    fig.tight_layout(); return fig


# ───────────────────────── 04 · clustering coefficient ──────────────────────────


# ───────────────────────── 04c · known-pair recovery ────────────────────────────
def fig04_recovery(recovery_results):
    ks = [1, 2, 3, 4, 5]
    fig, ax = plt.subplots(figsize=(7.4, 4.7))
    for kg in _present(recovery_results.keys()):
        r = recovery_results.get(kg, {}).get('recovery', {})
        ys = [r.get(k, np.nan) * 100 for k in ks]
        ax.plot(ks, ys, marker='o', ms=5, lw=1.6, color=kgc(kg), mfc='white', mec=kgc(kg),
                mew=1.4, label=NICE.get(kg, kg))
    ax.set_xticks(ks); ax.set_xlabel('hop distance (k)')
    ax.set_ylabel('% of known pairs recovered within k hops'); ax.set_ylim(0, 104)
    ax.legend(loc='lower right', ncol=2)
    _title(ax, 'Known drug–disease pair recovery by hop distance')
    fig.tight_layout(); return fig


# ───────────────────────── 04 · differential resilience (dropout curves) ────────


# ───────────────────────── 04 · differential resilience (dropout curves) ────────
def fig04_resilience(diversity_results):
    from matplotlib.lines import Line2D
    kgs = _present(diversity_results.keys())
    ncol = len(kgs)
    fig, axes = plt.subplots(1, ncol, figsize=(2.7 * ncol + 0.6, 4.2), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, kg in zip(axes, kgs):
        r = diversity_results[kg]
        x = [v * 100 for v in r['drop_rates']]
        gc, rc = r['gold_curve'], r['rand_curve']; c = kgc(kg)
        ax.plot(x, gc, '-o', color=c, lw=2, ms=5, label='gold standard')
        ax.plot(x, rc, '--s', color='#9A9A9A', lw=1.8, ms=5, label='random pairs')
        ax.fill_between(x, rc, gc, where=[g >= rr for g, rr in zip(gc, rc)], alpha=0.15, color=c)
        ax.text(0.96, 0.05, f"DR {r['dr_score']:.2f}", transform=ax.transAxes, ha='right', va='bottom',
                fontsize=8.5, fontweight='bold', color=TEXT_COLOR,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ddd'))
        ax.set_title(NICE.get(kg, kg), fontsize=10.5, fontweight='bold', color=TEXT_COLOR)
        ax.set_xlabel('edge dropout (%)', fontsize=8.5); ax.set_ylim(0, 1.05); ax.set_xlim(-2, max(x) + 2)
    axes[0].set_ylabel('recovery fraction', fontsize=9)
    axes[0].legend(handles=[Line2D([0], [0], color='#555', lw=2, marker='o', ms=5, label='gold standard'),
                            Line2D([0], [0], color='#9A9A9A', lw=1.8, ls='--', marker='s', ms=5, label='random pairs')],
                   loc='lower left', fontsize=8)
    fig.suptitle('Differential resilience to edge dropout: gold vs random pairs', fontsize=12,
                 fontweight='bold', x=0.012, ha='left', color=TEXT_COLOR, y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.94]); return fig


# ───────────────────────── 05 · ranking stability under dropout ─────────────────
