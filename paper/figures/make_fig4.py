"""Fig. 4 redrawn from results tables in the house style.

a  slope chart: link-prediction AUROC (type-constrained, Adamic-Adar) vs
   neighbourhood retrieval Recall@K at K = 0.5% of the pool (the scored variant)
b  graph-reachability index (Supp Table S3) vs link-prediction AUROC
c  AUROC on all diseases vs the sparse quartile, with % change and n positives
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])
sys.path.insert(0, str(ROOT / 'src'))
from biokgsuite.viz.style import (setup_style, clean_ax, panel_label, KG_ORDER, KG_NAMES,
                                  KG_PALETTE, NEUTRAL, FAINT, MUTED, TEXT_COLOR, fig_size)

setup_style()
T = ROOT / 'results' / 'tables'
lp = pd.read_csv(T / '06_task_performance' / 'link_prediction.csv')
lp = lp[(lp['heuristic'] == 'Adamic-Adar') & (lp['neg_strategy'] == 'type-constrained')]

auroc = dict(zip(lp['kg'], lp['auroc']))
nr = pd.read_csv(T / '06_task_performance' / 'neighbourhood_retrieval.csv')
recall = dict(zip(nr['KG'], nr['recall@frac']))
kfrac = dict(zip(nr['KG'], nr['k_norm'])) if 'k_norm' in nr.columns else {}
topo = pd.read_csv(T / '04_topology' / 'dimension_scores.csv')
reach = dict(zip(topo['KG'], topo['reachability']))
sp = pd.read_csv(T / '07_generalization' / 'sparse_tiers.csv')
q1 = sp[sp['tier'].str.startswith('Q1')]
sparse = dict(zip(q1['KG'], q1['auroc']))
npos = dict(zip(q1['KG'], q1['n_pos'].astype(int)))

kgs = KG_ORDER
fig, axes = plt.subplots(1, 3, figsize=(fig_size(2)[0], 2.45),
                         gridspec_kw={'width_ratios': [1.0, 1.0, 1.15], 'wspace': 0.75})

# ── a: slope chart ──────────────────────────────────────────────────────────
ax = axes[0]
rank_lp = {k: r for r, k in enumerate(sorted(kgs, key=lambda k: -auroc[k]))}
rank_nr = {k: r for r, k in enumerate(sorted(kgs, key=lambda k: -recall[k]))}
n = len(kgs)
for k in kgs:
    y0, y1 = n - rank_lp[k], n - rank_nr[k]
    hold = rank_lp[k] == rank_nr[k]
    ax.plot([0, 1], [y0, y1], color=KG_PALETTE[k], lw=1.0,
            ls=(0, (1.2, 1.2)) if hold else '-', zorder=2)
    ax.plot(0, y0, 'o', color=KG_PALETTE[k], ms=4, zorder=3)
    ax.plot(1, y1, 'o', color=KG_PALETTE[k], ms=4, zorder=3)
    ax.text(-0.08, y0, f'{KG_NAMES[k]}  {auroc[k]:.3f}', ha='right', va='center', fontsize=6, color=TEXT_COLOR)
    ax.text(1.08, y1, f'{recall[k]:.3f}', ha='left', va='center', fontsize=6, color=TEXT_COLOR)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(0.4, n + 0.6)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Link prediction\n(AUROC)', 'Gene retrieval\n(Recall@K, K = 0.5%)'], fontsize=6)
ax.set_yticks([])
for s in ('left', 'top', 'right'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', length=0)
ax.text(-0.08, n + 0.55, 'rank', ha='right', va='bottom', fontsize=6, color=MUTED)
panel_label(ax, 'a', x=-0.55)

# ── b: reachability vs AUROC ────────────────────────────────────────────────
ax = axes[1]
ax.axvspan(0.5, 1.02, color='#F2F2F2', zorder=0, lw=0)
for k in kgs:
    ax.plot(reach[k], auroc[k], 'o', color=KG_PALETTE[k], ms=5, zorder=3)
offsets = {'primekg': (4, 3), 'hetionet': (4, -7), 'drkg': (-4, -8), 'openbilink': (5, -2),
           'biokg': (-4, 4), 'matrix': (-4, -8)}
ha = {'primekg': 'left', 'hetionet': 'left', 'drkg': 'right', 'openbilink': 'left', 'biokg': 'right', 'matrix': 'right'}
for k in kgs:
    ax.annotate(KG_NAMES[k], (reach[k], auroc[k]), xytext=offsets[k], textcoords='offset points',
                fontsize=6, color=TEXT_COLOR, ha=ha[k], va='center')
ax.text(0.52, 0.985, 'no ordering\nabove 0.5', fontsize=6, color=MUTED, va='top', ha='left')
ax.set_xlim(0.38, 1.02)
ax.set_ylim(0.70, 1.01)
ax.set_xticks([0.4, 0.6, 0.8, 1.0])
ax.set_yticks([0.7, 0.8, 0.9, 1.0])
clean_ax(ax, xlabel='Graph reachability (index)', ylabel='Link-prediction AUROC')
panel_label(ax, 'b', x=-0.32)

# ── c: dumbbell all vs sparse ───────────────────────────────────────────────
ax = axes[2]
order = sorted(kgs, key=lambda k: -auroc[k])
for i, k in enumerate(order):
    y = len(order) - i
    a, s = auroc[k], sparse[k]
    ax.plot([s, a], [y, y], color=KG_PALETTE[k], lw=1.2, zorder=2)
    ax.plot(a, y, 'o', color=KG_PALETTE[k], ms=5, zorder=3)
    ax.plot(s, y, 'o', mfc='white', mec=KG_PALETTE[k], mew=1.0, ms=5, zorder=3)
    ax.text(0.47, y + 0.02, KG_NAMES[k], ha='right', va='bottom', fontsize=6, color=TEXT_COLOR)
    ax.text(0.47, y - 0.06, f'n = {npos[k]}', ha='right', va='top', fontsize=5.5, color=MUTED)
    ax.text(1.03, y, f'{100 * (s - a) / a:+.0f}%', ha='left', va='center', fontsize=6, color=TEXT_COLOR)
ax.axvline(0.5, color=FAINT, ls='--', lw=0.6, zorder=1)
ax.text(0.5, len(order) + 0.75, 'chance', fontsize=6, color=MUTED, ha='center', va='bottom')
ax.set_xlim(0.45, 1.02)
ax.set_ylim(0.4, len(order) + 0.9)
ax.set_xticks([0.5, 0.75, 1.0])
ax.set_yticks([])
clean_ax(ax, xlabel='Link-prediction AUROC')
ax.spines['left'].set_visible(False)
ax.plot([], [], 'o', color=TEXT_COLOR, ms=4, label='All diseases')
ax.plot([], [], 'o', mfc='white', mec=TEXT_COLOR, ms=4, label='Sparse quartile (Q1)')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=2, fontsize=6, frameon=False)
panel_label(ax, 'c', x=-0.42)

OUT.mkdir(parents=True, exist_ok=True)
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'Fig4_task_performance.{ext}', dpi=600, bbox_inches='tight', facecolor='white')
print('AUROC', {k: round(auroc[k], 3) for k in kgs})
print('Recall@K', {k: round(recall[k], 3) for k in kgs})
print('reach', {k: round(reach[k], 3) for k in kgs})
print('sparse', {k: (round(sparse[k], 3), npos[k]) for k in kgs})
print('hold rank:', [k for k in kgs if rank_lp[k] == rank_nr[k]])
