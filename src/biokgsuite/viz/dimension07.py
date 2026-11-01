"""Figures for dimension 07, generalisation.

One module per dimension, under the package rather than beside the notebook, so
there is a single import root and no `sys.path` manipulation. Each function
returns its figure rather than saving it; saving belongs to
`biokgsuite.io.Outputs`, so every notebook writes through one path and reports
what it produced.

Every function takes the graphs it is given (``kgs`` or the keys of a results
dict) and draws them in ``KG_ORDER``; the ``kg_colors`` arguments are kept for
the notebooks' call sites but a graph is always drawn in ``KG_PALETTE``.
"""
import textwrap

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from biokgsuite.viz.style import (setup_style, clean_ax, facet_title, fig_size,
                                  FAINT, HEATMAP_CMAP, MUTED, NEUTRAL, TEXT_COLOR)
from biokgsuite.viz.primitives import (kgc, _dark_cell, _present, KGS, NICE, facet_grid,
                                       panel_letters, ref_line)

setup_style()

BAR_W = 0.6          # of the slot, shared by every grouped bar chart in the paper
HEAT_ROW_IN = 0.28   # one heatmap row; fixed so two and six graphs keep the same cell size


def _ordered(names):
    """``names`` in KG_ORDER; anything that is not one of the six keeps its place after them."""
    names = list(names)
    return _present(names) + [n for n in names if n not in KGS]


def _empty(msg, figsize=None):
    fig, ax = plt.subplots(figsize=figsize or fig_size(1, 1), constrained_layout=True)
    ax.text(0.5, 0.5, msg, ha='center', va='center', transform=ax.transAxes,
            fontsize=7, color=MUTED)
    ax.set_axis_off()
    return fig


def _grouped_bars(ax, groups, kgs, values, show_values=False, fmt='{:.2f}', base=0.0):
    """Bars for ``kgs`` inside each of ``groups``; the group takes 0.6 of its slot.

    ``values[kg]`` is aligned with ``groups``. A NaN is an empty slot with a
    muted 'n/a' at the base line, so it cannot be read as a measured zero.
    """
    n = max(1, len(kgs))
    w = BAR_W / n
    x = np.arange(len(groups))
    for i, kg in enumerate(kgs):
        xs = x + (i - (n - 1) / 2) * w
        vals = np.asarray(values[kg], dtype=float)
        ax.bar(xs, np.where(np.isnan(vals), base, vals) - base, w, bottom=base,
               color=kgc(kg), label=NICE.get(kg, kg), zorder=2)
        for xi, v in zip(xs, vals):
            if np.isnan(v):
                ax.text(xi, base + 0.005 * (ax.get_ylim()[1] - base), 'n/a', ha='center',
                        va='bottom', fontsize=6, color=MUTED, style='italic', rotation=90)
            elif show_values:
                ax.text(xi, v + 0.01 * (ax.get_ylim()[1] - base), fmt.format(v),
                        ha='center', va='bottom', fontsize=6, color=TEXT_COLOR)
    ax.set_xticks(x)
    ax.set_xlim(-0.5, len(groups) - 0.5)


def _kg_legend(fig, kgs, ncol=None):
    """One legend for the figure, under the panels, graph colour to graph name."""
    handles = [Patch(color=kgc(k), label=NICE.get(k, k)) for k in kgs]
    fig.legend(handles=handles, loc='outside lower center', ncol=ncol or min(6, len(kgs)))


# ───────────────────────── 07a · sparse tiers ───────────────────────────────────
def fig07_sparse_tiers(tier_results, tiers, kgs, kg_colors):
    """AUROC and AUPRC by disease-connectivity tier, Q1 (sparse) highlighted."""
    shown = _ordered(n for n in kgs if tier_results.get(n))
    if not shown:
        return _empty('no tier results')
    labels = ['Q1\n(sparse)', 'Q2\n(under-studied)', 'Q3\n(moderate)', 'Q4\n(well-studied)']
    labels = labels[:len(tiers)] + [str(t) for t in tiers[len(labels):]]
    metrics = [('auroc', 'AUROC', 0.45), ('auprc', 'AUPRC', 0.0)]
    fig, axes = plt.subplots(1, 2, figsize=fig_size(2, 1), constrained_layout=True)
    for ax, (metric, label, ymin) in zip(axes, metrics):
        values = {n: [tier_results[n].get(t, {}).get(metric, np.nan) for t in tiers] for n in shown}
        ax.set_ylim(ymin, 1.0)
        # The sparse tier is the one the section is about; a faint band marks
        # it without competing with the bars.
        ax.axvspan(-0.5, 0.5, color=FAINT, alpha=0.25, lw=0, zorder=0)
        _grouped_bars(ax, tiers, shown, values, base=ymin)
        # AUROC's chance floor drawn, so a 0.52 bar is not read on the same
        # scale of success as a 0.95 one when the axis starts at 0.45.
        if metric == 'auroc':
            ref_line(ax, 0.5, 'chance')
        ax.set_xticklabels(labels)
        clean_ax(ax, ylabel=label)
        facet_title(ax, label)
    if len(shown) > 1:
        _kg_legend(fig, shown)
    panel_letters(fig, axes)
    return fig


# ───────────────────────── 07b · domain heatmap ─────────────────────────────────
def fig07_domain_heatmap(area_results, all_domains, kgs):
    """Per-domain AUROC, one row per graph. Unevaluable cells are grey n/a."""
    shown = _ordered(n for n in kgs if area_results.get(n))
    if not shown:
        return _empty('no domain results')
    areas = sorted(a for a in all_domains
                   if any(not np.isnan(area_results[n].get(a, {}).get('auroc', np.nan))
                          for n in shown))
    if not areas:
        return _empty('no domain cleared the minimum pair count')
    n_areas, n_kgs = len(areas), len(shown)
    vmin, vmax = 0.5, 1.0
    # Up to four domain names fit horizontally in a single column; more need
    # the double column. Long names wrap to two lines rather than rotate.
    width = fig_size(1)[0] if n_areas <= 4 else fig_size(2)[0]
    fig, ax = plt.subplots(figsize=(width, HEAT_ROW_IN * n_kgs + 0.6), constrained_layout=True)
    heat = np.array([[area_results[n].get(a, {}).get('auroc', np.nan)
                      for a in areas] for n in shown])
    valid = ~np.isnan(heat)
    im = ax.pcolormesh(np.where(valid, heat, np.nan), cmap=HEATMAP_CMAP, vmin=vmin, vmax=vmax,
                       edgecolors='none', linewidth=0, antialiased=False)
    ax.pcolormesh(np.where(~valid, 1, np.nan),
                  cmap=mcolors.ListedColormap(['#F2F2F2']), vmin=0, vmax=1,
                  edgecolors='none', linewidth=0, antialiased=False)
    ax.set_xticks([i + 0.5 for i in range(n_areas)])
    ax.set_xticklabels(['\n'.join(textwrap.wrap(str(a), 12, break_long_words=False)) for a in areas],
                       fontsize=6, color=TEXT_COLOR)
    ax.set_yticks([i + 0.5 for i in range(n_kgs)])
    ax.set_yticklabels([NICE.get(n, n) for n in shown], fontsize=6, color=TEXT_COLOR)
    ax.tick_params(length=0, which='both', pad=3)
    ax.grid(False, which='both')
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(n_kgs):
        for j in range(n_areas):
            v = heat[i, j]
            if valid[i, j]:
                frac = (v - vmin) / (vmax - vmin)
                ax.text(j + 0.5, i + 0.5, f'{v:.2f}', ha='center', va='center', fontsize=6,
                        color='white' if _dark_cell(HEATMAP_CMAP, frac) else TEXT_COLOR)
            else:
                ax.text(j + 0.5, i + 0.5, 'n/a', ha='center', va='center',
                        fontsize=6, color=MUTED, style='italic')
    # The bar lives in an inset so it is exactly as tall as the cells.
    cax = ax.inset_axes([1.02, 0, 0.03 if n_areas > 4 else 0.05, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=6, length=2, width=0.5, colors=TEXT_COLOR)
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(TEXT_COLOR)
    cbar.set_label('AUROC', fontsize=6, color=TEXT_COLOR, labelpad=3)
    return fig


# ───────────────────────── 07c · domain CV ──────────────────────────────────────
def fig07_domain_cv(area_cv, kgs, kg_colors, threshold=10.0):
    """Coefficient of variation of per-domain AUROC, one bar per graph.

    The message is which graphs clear the threshold, so the bars carry no
    values; the table beside the figure has them.
    """
    cv_kgs = _ordered(n for n in kgs if not np.isnan(area_cv.get(n, np.nan)))
    if not cv_kgs:
        return _empty('no CV values')
    vals = np.array([area_cv[n] for n in cv_kgs], dtype=float)
    xs = np.arange(len(cv_kgs))
    # Six graph names need the double column to stay horizontal at 6 pt.
    fig, ax = plt.subplots(figsize=fig_size(1 if len(cv_kgs) <= 4 else 2, 1),
                           constrained_layout=True)
    ax.bar(xs, vals, width=BAR_W, color=[kgc(n) for n in cv_kgs], zorder=2)
    ax.set_ylim(0, max(vals.max(), threshold) * 1.15)
    ax.set_xlim(-0.6, len(cv_kgs) - 0.4)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ref_line(ax, threshold, f'{threshold:.0f}% threshold')
    ax.set_xticks(xs)
    ax.set_xticklabels([NICE.get(n, n) for n in cv_kgs])
    clean_ax(ax, ylabel='CV of per-domain AUROC (%)')
    return fig


# ───────────────────────── 07e · prospective Hits@K ─────────────────────────────
def fig07_prospective_hits(temporal_metrics, kgs, kg_colors):
    """Hits@K on post-cutoff approvals, one group per K and one bar per graph."""
    keys = ['hits@10', 'hits@50', 'hits@100']
    labels = ['Hits@10', 'Hits@50', 'Hits@100']
    shown = _ordered(n for n in kgs if temporal_metrics.get(n, {}).get('hits@10') is not None)
    if not shown:
        return _empty('no graph cleared the minimum pair count')
    values = {n: [temporal_metrics.get(n, {}).get(k, np.nan) for k in keys] for n in shown}
    fig, ax = plt.subplots(figsize=fig_size(1, 1), constrained_layout=True)
    ax.set_ylim(0, 1.0)
    # The fraction is the message and six bars is the ceiling for printing it.
    _grouped_bars(ax, keys, shown, values, show_values=len(shown) * len(keys) <= 6,
                  fmt='{:.0%}')
    ax.set_xticklabels(labels)
    clean_ax(ax, ylabel='Fraction of new pairs')
    if len(shown) > 1:
        _kg_legend(fig, shown, ncol=3)
    return fig


# ───────────────────────── 07f · rank distribution ──────────────────────────────
def fig07_rank_distribution(pct_data, kg_colors, seed=42):
    """Rank-percentile distribution of new pairs, one row per graph.

    The vertical jitter is cosmetic but it was drawn from the GLOBAL numpy RNG,
    so the figure differed on every run and could not be reproduced from the
    same results. It now draws from its own seeded generator.

    Each row is a strip of pairs with the median as a short vertical tick; the
    median's value is written past the row end, so no legend is needed.
    """
    if not pct_data:
        return _empty('no rank data')
    rng = np.random.RandomState(seed)
    names = _ordered(pct_data.keys())
    fig, ax = plt.subplots(figsize=(fig_size(1)[0], 0.3 * len(names) + 0.55),
                           constrained_layout=True)
    all_pct = []
    for i, name in enumerate(names):
        pct = np.asarray(pct_data[name], dtype=float)
        if pct.size == 0:
            ax.annotate('no pairs', xy=(1, i), xycoords=('axes fraction', 'data'),
                        xytext=(3, 0), textcoords='offset points', ha='left', va='center',
                        fontsize=6, color=MUTED, style='italic')
            continue
        all_pct.extend(pct.tolist())
        med = float(np.median(pct))
        ax.scatter(pct, rng.normal(i, 0.08, size=len(pct)), color=kgc(name), alpha=0.6,
                   s=8, edgecolors='none', zorder=2)
        ax.plot([med, med], [i - 0.3, i + 0.3], color=TEXT_COLOR, lw=1, zorder=3)
        ax.annotate(f'median {med:.1f}', xy=(1, i), xycoords=('axes fraction', 'data'),
                    xytext=(3, 0), textcoords='offset points', ha='left', va='center',
                    fontsize=6, color=MUTED)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([NICE.get(n, n) for n in names])
    ax.set_ylim(len(names) - 0.5, -0.5)
    # max() of an empty sequence raised here when every graph had no evaluable
    # pair, which is exactly the case this figure is drawn to show.
    hi = max(all_pct) if all_pct else 100.0
    ax.set_xlim(-0.5, hi + 0.04 * max(hi, 1))
    ax.tick_params(axis='y', length=0)
    clean_ax(ax, xlabel='Rank percentile (%)')
    return fig


# ───────────────────────── 07g · network proximity ──────────────────────────────
def fig07_proximity(all_proximity):
    """Observed vs expected network distance for true drug-disease pairs.

    Points below the diagonal are closer than random expectation. Colour marks
    the pairs whose permutation test cleared significance; grey marks the rest,
    so the panel shows the whole distribution rather than only the wins.
    """
    kgs = [k for k in _present(all_proximity.keys()) if all_proximity.get(k)]
    if not kgs:
        return _empty('no proximity results')

    fig, axes, ncol = facet_grid(len(kgs))
    for idx, (ax, kg) in enumerate(zip(axes, kgs)):
        ok = [r for r in all_proximity[kg] if 'obs_dist' in r]
        rx = np.array([r['random_mean'] for r in ok], dtype=float)
        oy = np.array([r['obs_dist'] for r in ok], dtype=float)
        sig = np.array([bool(r.get('significant')) for r in ok], dtype=bool)
        lim = max(max(rx.tolist() + oy.tolist(), default=1), 1) * 1.05
        ax.plot([0, lim], [0, lim], ls='--', color=FAINT, lw=0.6, zorder=1)
        ax.scatter(rx[~sig], oy[~sig], s=8, color=FAINT, edgecolor='white', lw=0.3, zorder=2)
        ax.scatter(rx[sig], oy[sig], s=10, color=kgc(kg), edgecolor='white', lw=0.3, zorder=3)
        frac = 100 * sig.sum() / len(ok) if ok else 0
        # n is the number of pairs that produced a distance at all, which is
        # smaller than the number attempted whenever an endpoint is missing or
        # unreachable. Shown so the percentage is read against its own base.
        # Top left is above the diagonal, where the significant pairs are not.
        ax.text(0.04, 0.96, f'{frac:.0f}% significant\n(n = {len(ok)})',
                transform=ax.transAxes, va='top', fontsize=6, color=MUTED)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect('equal')
        # Equal-aspect panels are narrow, so the axis label is short; the
        # caption says the expectation is the mean over random diseases.
        clean_ax(ax, xlabel='Expected distance' if idx >= len(kgs) - ncol else '',
                 ylabel='Observed distance' if idx % ncol == 0 else '')
        facet_title(ax, NICE.get(kg, kg))

    # The significant series is a different colour in every facet, so its
    # legend swatch is neutral ink; the facet title already says which graph.
    handles = [Line2D([0], [0], ls='', marker='o', ms=3.5, color=TEXT_COLOR,
                      label='Closer than random (permutation test)'),
               Line2D([0], [0], ls='', marker='o', ms=3.5, color=FAINT, label='Not significant')]
    fig.legend(handles=handles, loc='outside lower center', ncol=2)
    if len(kgs) > 1:
        panel_letters(fig, axes)
    return fig


# ───────────────────────── 07d · prospective AUROC ──────────────────────────────
def fig07_prospective_auroc(baseline_all, baseline_same, new_pairs):
    """Prospective AUROC on post-cutoff approvals, against two baselines.

    Each group on the x axis is one graph. The two baselines are greys (they
    are references, not graphs); the new-pairs bar wears the graph's colour.
    """
    kgs = [k for k in _present(new_pairs.keys())
           if not np.isnan(new_pairs.get(k, {}).get('auroc', np.nan))]
    if not kgs:
        return _empty('no prospective results')

    series = [('Baseline (all)', baseline_all, lambda k: NEUTRAL),
              ('Baseline (same)', baseline_same, lambda k: FAINT),
              ('New pairs', new_pairs, kgc)]
    x = np.arange(len(kgs))
    w = BAR_W / len(series)
    ymin = 0.45
    fig, ax = plt.subplots(figsize=fig_size(1 if len(kgs) <= 3 else 2, 1), constrained_layout=True)
    for i, (label, src, colour) in enumerate(series):
        vals = np.array([src.get(k, {}).get('auroc', np.nan) for k in kgs], dtype=float)
        xs = x + (i - 1) * w
        ax.bar(xs, np.where(np.isnan(vals), ymin, vals) - ymin, w, bottom=ymin,
               color=[colour(k) for k in kgs], zorder=2)
        for xi, v in zip(xs, vals):
            if np.isnan(v):
                # Just above the chance line, so the two do not cross.
                ax.text(xi, 0.515, 'n/a', ha='center', va='bottom', fontsize=6,
                        color=MUTED, style='italic', rotation=90)

    # Chance, drawn. The y axis starts at 0.45, so without this line a bar at
    # 0.52 and a bar at 0.95 look like points on the same scale of success.
    ref_line(ax, 0.5, 'chance')
    ax.set_xticks(x)
    ax.set_xlim(-0.5, len(kgs) - 0.5)
    ax.set_xticklabels([NICE.get(k, k) for k in kgs])
    ax.set_ylim(ymin, 1.0)
    clean_ax(ax, ylabel='AUROC')
    handles = [Patch(color=NEUTRAL, label='Baseline (all)'),
               Patch(color=FAINT, label='Baseline (same)'),
               Patch(color=TEXT_COLOR, label='New pairs (graph colour)')]
    fig.legend(handles=handles, loc='outside lower center', ncol=3)
    return fig
