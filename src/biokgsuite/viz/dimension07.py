"""Figures for dimension 07, generalisation.

One module per dimension, under the package rather than beside the notebook, so
there is a single import root and no `sys.path` manipulation. Every dimension
folder used to define a module named `figures`, and Python's module cache hands
whichever was imported first to every notebook that runs afterwards in the same
kernel, so `from figures import ...` in 07 could silently return 03's.

Each function RETURNS its figure. Saving belongs to `biokgsuite.io.Outputs`, so
every notebook writes through one path and reports what it produced. The former
`figures.py` called `_save()` itself, which is why these two figures never
appeared in the results manifest.

The file it replaces also carried a duplicated section header and a trailing
`03 · source-database diversity` header with no function under it, left over
from the split of the old shared `eval_figs.py`.
"""
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from biokgsuite.viz.style import (setup_style, clean_ax, panel_label,
                                  KG_PALETTE, TEXT_COLOR, ALERT_RED,
                                  HEATMAP_CMAP)
from biokgsuite.viz.primitives import kgc, _present, NICE

setup_style()


def _style_left_spine(ax):
    """Left spine only. Was defined inline in the notebook, so the figures could
    not be redrawn without executing the analysis cells above them."""
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.spines['left'].set_visible(True)


def _empty(msg, figsize=(4, 3)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.text(0.5, 0.5, msg, ha='center', va='center', transform=ax.transAxes,
            color=TEXT_COLOR)
    ax.set_axis_off()
    return fig


def fig07_sparse_tiers(tier_results, tiers, kgs, kg_colors):
    """AUROC and AUPRC by disease-connectivity tier, Q1 (sparse) highlighted."""
    shown = [n for n in kgs if tier_results.get(n)]
    if not shown:
        return _empty('no tier results')
    labels = ['Q1\n(Sparse)', 'Q2\n(Under-studied)', 'Q3\n(Moderate)',
              'Q4\n(Well-studied)']
    metrics = [('auroc', 'AUROC', 0.45), ('auprc', 'AUPRC', 0.0)]
    x = np.arange(len(tiers))
    w = 0.15
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, (metric, label, ymin), letter in zip(axes, metrics, 'ab'):
        for i, name in enumerate(shown):
            vals = [tier_results[name].get(t, {}).get(metric, np.nan) for t in tiers]
            bars = ax.bar(x + i * w, vals, w, label=NICE.get(name, name.upper()),
                          color=kg_colors[name], edgecolor='none', linewidth=0)
            for bar, v in zip(bars, vals):
                if not np.isnan(v):
                    ax.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_height() + 0.01, f'{v:.2f}', ha='center',
                            va='bottom', fontsize=7, color=TEXT_COLOR, rotation=90)
        ax.axvspan(-0.4, 0.4 + w * (len(shown) - 1), alpha=0.06,
                   color=ALERT_RED, zorder=0)
        # AUROC's chance floor drawn, so a 0.52 bar is not read on the same
        # scale of success as a 0.95 one when the axis starts at 0.45.
        if metric == 'auroc':
            ax.axhline(0.5, ls=':', lw=1, color='#999', zorder=1)
        ax.set_xticks(x + w * (len(shown) - 1) / 2)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylim(ymin, 1.15)
        ax.set_ylabel(label, fontsize=10)
        ax.legend(fontsize=9, frameon=True, facecolor='white', edgecolor='none',
                  framealpha=0.9, loc='upper left')
        _style_left_spine(ax)
        clean_ax(ax, title=f'{label} by Disease Connectivity Tier')
        panel_label(ax, letter)
    fig.suptitle('Data-Sparse Generalisation', fontsize=12, fontweight='bold',
                 y=1.02, color=TEXT_COLOR)
    fig.tight_layout()
    return fig


def fig07_domain_heatmap(area_results, all_domains, kgs):
    """Per-domain AUROC, one row per graph. Unevaluable cells are grey N/A."""
    shown = [n for n in kgs if area_results.get(n)]
    if not shown:
        return _empty('no domain results')
    areas = sorted(a for a in all_domains
                   if any(not np.isnan(area_results[n].get(a, {}).get('auroc', np.nan))
                          for n in shown))
    if not areas:
        return _empty('no domain cleared the minimum pair count')
    n_areas, n_kgs = len(areas), len(shown)
    fig, ax = plt.subplots(figsize=(max(6.0, n_areas * 0.72),
                                    max(2.4, n_kgs * 0.72 + 0.8)))
    heat = np.array([[area_results[n].get(a, {}).get('auroc', np.nan)
                      for a in areas] for n in shown])
    im = ax.pcolormesh(heat, cmap=HEATMAP_CMAP, vmin=0.5, vmax=1.0,
                       edgecolors='none', linewidth=0, antialiased=False)
    ax.pcolormesh(np.where(np.isnan(heat), 1, np.nan),
                  cmap=mcolors.ListedColormap(['#e8e8e8']), vmin=0, vmax=1,
                  edgecolors='none', linewidth=0, antialiased=False)
    ax.set_xticks([i + 0.5 for i in range(n_areas)])
    ax.set_xticklabels(areas, rotation=40, ha='right', fontsize=8.5)
    ax.set_yticks([i + 0.5 for i in range(n_kgs)])
    ax.set_yticklabels([NICE.get(n, n.upper()) for n in shown], fontsize=9)
    ax.tick_params(length=0, which='both')
    ax.grid(False, which='both')
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(n_kgs):
        for j in range(n_areas):
            v = heat[i, j]
            if not np.isnan(v):
                ax.text(j + 0.5, i + 0.5, f'{v:.2f}', ha='center', va='center',
                        fontsize=8, fontweight='bold',
                        color='white' if v > 0.72 else TEXT_COLOR)
            else:
                ax.text(j + 0.5, i + 0.5, 'N/A', ha='center', va='center',
                        fontsize=7, color='#aaaaaa', style='italic')
    cbar = fig.colorbar(im, ax=ax, shrink=0.80, pad=0.02, label='AUROC')
    cbar.ax.tick_params(labelsize=7.5)
    ax.set_title('AUROC by domain', fontsize=10, fontweight='bold',
                 color=TEXT_COLOR, pad=10)
    panel_label(ax, 'a')
    fig.tight_layout()
    return fig


def fig07_domain_cv(area_cv, kgs, kg_colors, threshold=10.0):
    """Coefficient of variation of per-domain AUROC, one bar per graph."""
    cv_kgs = [n for n in kgs if not np.isnan(area_cv.get(n, np.nan))]
    if not cv_kgs:
        return _empty('no CV values')
    vals = [area_cv[n] for n in cv_kgs]
    fig, ax = plt.subplots(figsize=(3.0, 3.4))
    bars = ax.bar(range(len(cv_kgs)), vals,
                  color=[kg_colors[n] for n in cv_kgs], edgecolor='white',
                  linewidth=0.5)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                f'{v:.1f}%', ha='center', fontsize=9, fontweight='bold',
                color=TEXT_COLOR)
    ax.axhline(threshold, color=ALERT_RED, ls='--', lw=1.2, alpha=0.85,
               label=f'{threshold:.0f}% threshold')
    ax.set_xticks(range(len(cv_kgs)))
    ax.set_xticklabels([NICE.get(n, n.upper()) for n in cv_kgs], fontsize=8.5)
    ax.set_ylim(0, max(vals + [threshold * 1.2]) * 1.28)
    ax.legend(fontsize=7.5, frameon=False)
    _style_left_spine(ax)
    clean_ax(ax, title='AUROC CV\n(domain equity)', ylabel='CV (%)')
    panel_label(ax, 'c')
    fig.tight_layout()
    return fig


def fig07_prospective_hits(temporal_metrics, kgs, kg_colors):
    """Hits@K on post-cutoff approvals."""
    keys = ['hits@10', 'hits@50', 'hits@100']
    labels = ['Hits@10', 'Hits@50', 'Hits@100']
    shown = [n for n in kgs if temporal_metrics.get(n, {}).get('hits@10') is not None]
    if not shown:
        return _empty('no graph cleared the minimum pair count')
    x = np.arange(len(keys))
    w = 0.30
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    for i, name in enumerate(shown):
        m = temporal_metrics.get(name, {})
        vals = [m.get(k, np.nan) for k in keys]
        # A missing value is drawn as a zero-height bar, which reads as a
        # measured zero. Labelled below so the two are distinguishable.
        bars = ax.bar(x + i * w, [0 if np.isnan(v) else v for v in vals], w,
                      label=NICE.get(name, name.upper()), color=kg_colors[name],
                      edgecolor='none', linewidth=0)
        for bar, v in zip(bars, vals):
            if np.isnan(v):
                ax.text(bar.get_x() + bar.get_width() / 2, 0.01, 'NA',
                        ha='center', va='bottom', fontsize=8, color='#999')
            elif v > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.01, f'{v:.0%}', ha='center',
                        va='bottom', fontsize=9, color=TEXT_COLOR)
    ax.set_xticks(x + w * (len(shown) - 1) / 2)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel('Fraction', fontsize=10)
    ax.legend(fontsize=9, frameon=True, facecolor='white', edgecolor='none',
              framealpha=0.9)
    _style_left_spine(ax)
    clean_ax(ax, title='Hits@K (new pairs)')
    fig.tight_layout()
    return fig


def fig07_rank_distribution(pct_data, kg_colors, seed=42):
    """Rank-percentile distribution of new pairs, one row per graph.

    The vertical jitter is cosmetic but it was drawn from the GLOBAL numpy RNG,
    so the figure differed on every run and could not be reproduced from the
    same results. It now draws from its own seeded generator.
    """
    if not pct_data:
        return _empty('no rank data')
    rng = np.random.RandomState(seed)
    fig, ax = plt.subplots(figsize=(5, 2.5))
    all_pct = []
    for i, (name, pct) in enumerate(pct_data.items()):
        pct = np.asarray(pct, dtype=float)
        if pct.size == 0:
            continue
        all_pct.extend(pct.tolist())
        med = float(np.median(pct))
        ax.scatter(pct, rng.normal(i, 0.08, size=len(pct)),
                   color=kg_colors[name], alpha=0.6, s=25, edgecolors='none',
                   label=f'{NICE.get(name, name.upper())} (median: {med:.1f}%ile)')
        ax.axvline(med, color=kg_colors[name], ls='--', lw=1.2, alpha=0.5)
    ax.set_yticks(range(len(pct_data)))
    ax.set_yticklabels([NICE.get(n, n.upper()) for n in pct_data], fontsize=10)
    ax.set_xlabel('Rank percentile', fontsize=10)
    # max() of an empty sequence raised here when every graph had no evaluable
    # pair, which is exactly the case this figure is drawn to show.
    ax.set_xlim(-0.2, (max(all_pct) if all_pct else 100.0) + 0.5)
    _style_left_spine(ax)
    clean_ax(ax, title='Rank distribution (new pairs)')
    ax.legend(fontsize=9, frameon=True, facecolor='white', edgecolor='none',
              framealpha=0.9)
    fig.tight_layout()
    return fig


def fig07_proximity(all_proximity):
    """Observed vs expected network distance for true drug-disease pairs.

    Points below the diagonal are closer than random expectation. Colour marks
    the pairs whose permutation test cleared significance; grey marks the rest,
    so the panel shows the whole distribution rather than only the wins.
    """
    kgs = [k for k in _present(all_proximity.keys()) if all_proximity.get(k)]
    if not kgs:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, 'no proximity results', ha='center', va='center',
                transform=ax.transAxes, color=TEXT_COLOR)
        ax.set_axis_off()
        return fig

    ncol = min(3, len(kgs))
    nrow = math.ceil(len(kgs) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.9 * ncol, 3.5 * nrow),
                             sharex=True, sharey=True)
    axes = np.array(axes).reshape(-1)

    for ax, kg in zip(axes, kgs):
        ok = [r for r in all_proximity[kg] if 'obs_dist' in r]
        rx = [r['random_mean'] for r in ok]
        oy = [r['obs_dist'] for r in ok]
        sig = [bool(r.get('significant')) for r in ok]
        lim = max(max(rx + oy, default=1), 1) * 1.05
        ax.plot([0, lim], [0, lim], ls='--', color='#bbb', lw=1, zorder=1)
        ax.scatter([rx[i] for i in range(len(ok)) if not sig[i]],
                   [oy[i] for i in range(len(ok)) if not sig[i]],
                   s=16, color='#CCCCCC', edgecolor='white', lw=0.3, zorder=2,
                   label='n.s.')
        ax.scatter([rx[i] for i in range(len(ok)) if sig[i]],
                   [oy[i] for i in range(len(ok)) if sig[i]],
                   s=20, color=kgc(kg), edgecolor='white', lw=0.3, zorder=3,
                   label='closer than random')
        frac = 100 * sum(sig) / len(ok) if ok else 0
        # n is the number of pairs that produced a distance at all, which is
        # smaller than the number attempted whenever an endpoint is missing or
        # unreachable. Shown so the percentage is read against its own base.
        ax.text(0.04, 0.96, f'{frac:.0f}% significant\n(n={len(ok)})',
                transform=ax.transAxes, va='top', fontsize=8, color=TEXT_COLOR)
        ax.set_title(NICE.get(kg, kg), fontsize=10.5, fontweight='bold',
                     loc='left', color=TEXT_COLOR, pad=6)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect('equal')

    for ax in axes[len(kgs):]:
        ax.set_visible(False)

    fig.supxlabel('expected distance (random mean)', fontsize=9.5)
    fig.supylabel('observed distance', fontsize=9.5)
    axes[0].legend(loc='lower right', fontsize=7.5)
    fig.suptitle('Network proximity of true drug-disease pairs vs random expectation',
                 fontsize=12, fontweight='bold', x=0.012, ha='left',
                 color=TEXT_COLOR, y=1.0)
    fig.text(0.012, 0.945,
             'points below the diagonal are closer than random expectation',
             fontsize=9, color='#777', ha='left')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


def fig07_prospective_auroc(baseline_all, baseline_same, new_pairs):
    """Prospective AUROC on post-cutoff approvals, against two baselines.

    Series colours come from the KG palette purely as three distinguishable
    hues; they do not denote a graph. Each group on the x axis is one graph.
    """
    kgs = [k for k in _present(new_pairs.keys())
           if not np.isnan(new_pairs.get(k, {}).get('auroc', np.nan))]
    if not kgs:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, 'no prospective results', ha='center', va='center',
                transform=ax.transAxes, color=TEXT_COLOR)
        ax.set_axis_off()
        return fig

    series = [('Baseline (all)', baseline_all, KG_PALETTE['drkg']),
              ('Baseline (same)', baseline_same, KG_PALETTE['primekg']),
              ('New pairs', new_pairs, KG_PALETTE['hetionet'])]
    x = np.arange(len(kgs))
    w = 0.26
    fig, ax = plt.subplots(figsize=(1.5 * len(kgs) + 2.2, 4.3))
    for i, (label, src, color) in enumerate(series):
        vals = [src.get(k, {}).get('auroc', np.nan) for k in kgs]
        ax.bar(x + (i - 1) * w, vals, w, label=label, color=color,
               edgecolor='white', linewidth=0.5, zorder=2)

    # Chance, drawn. The y axis starts at 0.45, so without this line a bar at
    # 0.52 and a bar at 0.95 look like points on the same scale of success.
    ax.axhline(0.5, ls=':', lw=1, color='#999', zorder=1)
    ax.text(len(kgs) - 0.5, 0.503, 'chance', fontsize=7.5, color='#777',
            ha='right', va='bottom')

    ax.set_xticks(x)
    ax.set_xticklabels([NICE.get(k, k) for k in kgs])
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel('AUROC')
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=8.5)
    ax.set_title('Prospective AUROC: new pairs vs baselines', fontsize=11.5,
                 fontweight='bold', loc='left', color=TEXT_COLOR, pad=30)
    fig.tight_layout()
    return fig
