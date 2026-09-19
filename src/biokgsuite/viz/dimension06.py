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
import matplotlib.colors as mcolors

from biokgsuite.viz.style import (clean_ax, facet_title, fig_size, FAINT, HEATMAP_CMAP,
                                  HEURISTIC_COLORS, MUTED, NEUTRAL, TEXT_COLOR)
from biokgsuite.viz.primitives import (kgc, _dark_cell, _present, NICE, facet_grid,
                                       panel_letters)

STRATEGY_LABEL = {'random': 'Random', 'type-constrained': 'Type-constrained',
                  'shared-target': 'Shared-target'}
STRATEGY_ORDER = ['random', 'type-constrained', 'shared-target']

BAR_W = 0.6          # of the slot, shared by every grouped bar chart in the paper
HEAT_ROW_IN = 0.28   # one heatmap row; fixed so two and six graphs keep the same cell size


def _heuristic_label(name):
    """Scorer names in sentence case and British spelling; Adamic-Adar is a
    proper name and stays."""
    if '-' in name:
        return name
    return name.capitalize().replace('neighbors', 'neighbours')


def _grouped_bars(ax, groups, kgs, values, show_values):
    """Bars for ``kgs`` inside each of ``groups``; the group takes 0.6 of its slot.

    ``values[kg]`` is a sequence aligned with ``groups``. A NaN is drawn as an
    empty slot with a muted 'n/a', so it cannot be read as a measured zero.
    """
    n = max(1, len(kgs))
    w = BAR_W / n
    x = np.arange(len(groups))
    for i, kg in enumerate(kgs):
        xs = x + (i - (n - 1) / 2) * w
        vals = np.asarray(values[kg], dtype=float)
        ax.bar(xs, np.where(np.isnan(vals), 0, vals), w, color=kgc(kg),
               label=NICE.get(kg, kg), zorder=2)
        for xi, v in zip(xs, vals):
            if np.isnan(v):
                ax.text(xi, 0.005, 'n/a', ha='center', va='bottom', fontsize=6,
                        color=MUTED, style='italic', rotation=90)
            elif show_values:
                ax.text(xi, v + 0.01 * ax.get_ylim()[1], f'{v:.2f}', ha='center',
                        va='bottom', fontsize=6, color=TEXT_COLOR)
    ax.set_xticks(x)
    ax.set_xlim(-0.5, len(groups) - 0.5)


def _fraction_ylim(values):
    """Head room above the tallest bar, never past 1 (these are fractions)."""
    top = np.nanmax([np.nanmax(np.asarray(v, dtype=float)) for v in values] + [0.0])
    return min(1.05, max(0.1, top * 1.3))


def _place_note(fig, ax, text, curves, y_floor):
    """A muted note in whichever lower corner the curves leave empty.

    Candidates are tried left to right along the bottom of the panel, just
    above ``y_floor`` (the prevalence line); the first whose box no curve point
    enters wins, otherwise the one with the fewest points inside. A fixed corner
    cannot work for every graph: a precision-recall curve dips at low recall
    when the first positives rank late, and falls at high recall always.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    pts = np.vstack([ax.transData.transform(np.column_stack([np.asarray(r, float), np.asarray(p, float)]))
                     for r, p in curves]) if curves else np.empty((0, 2))
    best, best_hits = None, None
    # The upper-right corner is tried first: a precision-recall curve falls
    # towards the lower right, so high recall at high precision is the region
    # it enters least. The bottom candidates follow, just above the prevalence
    # line, for the rare panel whose curve stays high all the way.
    candidates = [((0.98, 0.98), 'right', 'top', ax.transAxes)]
    candidates += [((x, y_floor), ha, 'bottom', ax.get_xaxis_transform())
                   for x, ha in ((0.02, 'left'), (0.5, 'center'), (0.98, 'right'))]
    for (x, y), ha, va, transform in candidates:
        t = ax.text(x, y, text, transform=transform, ha=ha, va=va,
                    fontsize=6, color=MUTED, linespacing=1.3)
        bb = t.get_window_extent(renderer)
        hits = int(np.sum((pts[:, 0] >= bb.x0) & (pts[:, 0] <= bb.x1) &
                          (pts[:, 1] >= bb.y0) & (pts[:, 1] <= bb.y1))) if len(pts) else 0
        if best is None or hits < best_hits:
            if best is not None:
                best.remove()
            best, best_hits = t, hits
        else:
            t.remove()
        if best_hits == 0:
            break


# ───────────────────────── 06a · link prediction, heatmap ───────────────────────
def fig06_link_heatmap(link_records, primary_heuristic):
    """AUROC and AUPRC per KG per negative-sampling strategy, primary heuristic.

    Two facets of one table (one per metric) on a shared colour scale. The
    scorer is named once, in a muted note under the panels, because the panels
    are about the sampling strategies, not the scorer.
    """
    df = pd.DataFrame(link_records)
    df = df[df['heuristic'] == primary_heuristic]
    strat = [s for s in STRATEGY_ORDER if s in set(df['strategy'])] or list(df['strategy'].unique())
    kgs = _present(df['kg'].unique())
    vmin, vmax = 0.45, 1.0

    fig, axes = plt.subplots(1, 2, figsize=(fig_size(2)[0], HEAT_ROW_IN * max(1, len(kgs)) + 0.65),
                             constrained_layout=True)
    im = None
    for ax, metric, label in zip(axes, ['auroc', 'auprc'], ['AUROC', 'AUPRC']):
        piv = (df.pivot_table(index='kg', columns='strategy', values=metric, aggfunc='mean')
                 .reindex(index=kgs, columns=strat))
        data = piv.values.astype(float)
        valid = ~np.isnan(data)
        im = ax.pcolormesh(np.where(valid, data, np.nan), cmap=HEATMAP_CMAP,
                           vmin=vmin, vmax=vmax, edgecolors='none', linewidth=0,
                           antialiased=False)
        # "not evaluated" is a flat grey cell, distinct from a low score.
        ax.pcolormesh(np.where(~valid, 1, np.nan),
                      cmap=mcolors.ListedColormap(['#F2F2F2']), vmin=0, vmax=1,
                      edgecolors='none', linewidth=0, antialiased=False)
        ax.set_xticks([x + 0.5 for x in range(len(strat))])
        ax.set_xticklabels([STRATEGY_LABEL.get(s, s) for s in strat], fontsize=6, color=TEXT_COLOR)
        ax.set_yticks([y + 0.5 for y in range(len(kgs))])
        ax.set_yticklabels([NICE.get(k, k) for k in kgs], fontsize=6, color=TEXT_COLOR)
        ax.tick_params(length=0, pad=3)
        ax.invert_yaxis()
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
        for i in range(len(kgs)):
            for j in range(len(strat)):
                v = data[i, j]
                if valid[i, j]:
                    frac = (v - vmin) / (vmax - vmin)
                    ax.text(j + 0.5, i + 0.5, f'{v:.3f}', ha='center', va='center',
                            fontsize=6,
                            color='white' if _dark_cell(HEATMAP_CMAP, frac) else TEXT_COLOR)
                else:
                    ax.text(j + 0.5, i + 0.5, 'n/a', ha='center', va='center',
                            fontsize=6, color=MUTED, style='italic')
        facet_title(ax, label)

    # One bar for both panels: they share the scale, so a second bar would only
    # invite the reader to look for a difference that is not there.
    cax = axes[1].inset_axes([1.04, 0, 0.04, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=6, length=2, width=0.5, colors=TEXT_COLOR)
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(TEXT_COLOR)
    fig.supxlabel(f'{primary_heuristic} scorer', fontsize=6, color=MUTED, x=0.0, ha='left')
    panel_letters(fig, axes)
    return fig


# ───────────────────────── 06b · link prediction, PR curves ─────────────────────
def fig06_link_pr(pr_curves, strategy_label='type-constrained'):
    """Precision-recall per KG, one line per heuristic, against prevalence.

    `pr_curves` maps kg -> {'baseline': float,
                            'curves': {heuristic: (recall, precision, auprc)}}
    built from the SAME scores the AUPRC was computed from, so a curve and its
    label cannot disagree.

    ``strategy_label`` names the negative-sampling strategy for the caption; it
    is not drawn.
    """
    kgs = _present(pr_curves.keys())
    fig, axes, ncol = facet_grid(max(1, len(kgs)))
    handles = {}

    notes = []
    for idx, (ax, name) in enumerate(zip(axes, kgs)):
        entry = pr_curves[name]
        curves = []
        for heuristic, (rec, prec, _auprc) in entry['curves'].items():
            (line,) = ax.plot(rec, prec, color=HEURISTIC_COLORS.get(heuristic, NEUTRAL),
                              lw=1, label=_heuristic_label(heuristic), zorder=3)
            handles.setdefault(heuristic, line)
            curves.append((rec, prec))
        base = entry['baseline']
        ax.axhline(base, color=FAINT, ls='--', lw=0.6, zorder=1)
        # Below the line at the right end is the one region a precision-recall
        # curve never enters, so the label cannot land on a curve.
        ax.text(1.0, base - 0.02, f'prevalence {base:.3f}', ha='right', va='top',
                fontsize=6, color=MUTED)
        # AUPRC and prevalence are PER PANEL and belong in the panel. The shared
        # legend used to carry the first KG's numbers while sitting under every
        # panel, so Hetionet's AUPRC read as if it described BioKG too.
        lines = [f'{_heuristic_label(h)} {a:.3f}' for h, (_, _, a) in entry['curves'].items()]
        notes.append((ax, 'AUPRC\n' + '\n'.join(lines), curves, base + 0.03))
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.05)
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        clean_ax(ax, xlabel='Recall' if idx >= len(kgs) - ncol else '',
                 ylabel='Precision' if idx % ncol == 0 else '')
        if idx % ncol:
            ax.tick_params(labelleft=False)
        facet_title(ax, NICE.get(name, name))
    # Placed after every panel is laid out, since placement measures the render.
    for ax, text, curves, floor in notes:
        _place_note(fig, ax, text, curves, floor)

    order = [h for h in HEURISTIC_COLORS if h in handles] + [h for h in handles if h not in HEURISTIC_COLORS]
    fig.legend([handles[h] for h in order], [_heuristic_label(h) for h in order],
               loc='outside lower center', ncol=len(order))
    if len(kgs) > 1:
        panel_letters(fig, axes)
    return fig


# ───────────────────────── 06c · neighbourhood retrieval ────────────────────────
def fig06_retrieval(target_results, k_range=(1, 5, 10, 20, 50, 100, 200, 500)):
    """Recall@K curves, and Recall@10/50/100 as grouped bars."""
    kgs = _present([n for n in target_results if 'ranks' in target_results[n]])
    k_range = list(k_range)
    fig, axes = plt.subplots(1, 2, figsize=fig_size(2, 1), constrained_layout=True)

    ax = axes[0]
    curves = {}
    for name in kgs:
        ranks = np.asarray(target_results[name]['ranks'])
        ys = [float(np.mean(ranks <= k)) for k in k_range]
        curves[name] = ys
        ax.plot(k_range, ys, '-o', color=kgc(name), lw=1, ms=3.5, mfc='white',
                mec=kgc(name), mew=0.75, label=NICE.get(name, name), zorder=3)
    ax.set_ylim(0, _fraction_ylim(curves.values()) if curves else 1.05)
    ax.set_xlim(-0.02 * max(k_range), max(k_range) * 1.05)
    ax.set_xticks([t for t in ax.get_xticks() if 0 <= t <= max(k_range)])
    clean_ax(ax, xlabel='K', ylabel='Recall@K')
    facet_title(ax, 'Recall against K')

    ax = axes[1]
    ks = ['recall@10', 'recall@50', 'recall@100']
    values = {name: [target_results[name][k] for k in ks] for name in kgs}
    ax.set_ylim(0, _fraction_ylim(values.values()) if values else 1.05)
    _grouped_bars(ax, ks, kgs, values, show_values=len(kgs) * len(ks) <= 6)
    ax.set_xticklabels(['Recall@10', 'Recall@50', 'Recall@100'])
    clean_ax(ax, ylabel='Recall')
    facet_title(ax, 'Recall at fixed K')

    # One legend serves both panels; below the axes it can never sit on a
    # curve, wherever the six graphs happen to land.
    if len(kgs) > 1:
        fig.legend(*axes[0].get_legend_handles_labels(), loc='outside lower center',
                   ncol=min(6, len(kgs)))
    panel_letters(fig, axes)
    return fig


# ───────────────────────── 06d · multi-hop reasoning ────────────────────────────
def fig06_multihop(multihop_results):
    """Hits@K for the intermediate target, and the rank CDF behind it."""
    kgs = _present([n for n, r in multihop_results.items()
                    if r.get('n_evaluated', 0) > 0])
    fig, axes = plt.subplots(1, 2, figsize=fig_size(2, 1), constrained_layout=True)

    ax = axes[0]
    ks = ['hits@10', 'hits@50', 'hits@100']
    values = {name: [multihop_results[name][k] for k in ks] for name in kgs}
    ax.set_ylim(0, _fraction_ylim(values.values()) if values else 1.05)
    _grouped_bars(ax, ks, kgs, values, show_values=False)
    ax.set_xticklabels(['Hits@10', 'Hits@50', 'Hits@100'])
    clean_ax(ax, ylabel='Fraction recovered')
    facet_title(ax, 'Intermediate target recovery')

    ax = axes[1]
    for name in kgs:
        r = np.sort(np.asarray(multihop_results[name]['ranks']))
        ax.plot(r, np.arange(1, len(r) + 1) / len(r), '-', color=kgc(name), lw=1,
                label=f"{NICE.get(name, name)}, MRR {multihop_results[name]['mrr']:.3f}",
                zorder=3)
    ax.set_xscale('log')
    ax.set_ylim(0, 1.05)
    # K = 100 is the threshold behind Hits@100 in panel a. The label sits at the
    # top end of the line, above the axes, where no curve can reach it; the
    # facet title is kept short so the two never meet.
    ax.axvline(100, color=FAINT, ls='--', lw=0.6, zorder=1)
    ax.annotate('K = 100', xy=(100, 1), xycoords=('data', 'axes fraction'),
                xytext=(0, 2), textcoords='offset points',
                ha='center', va='bottom', fontsize=6, color=MUTED)
    clean_ax(ax, xlabel='Rank of known target', ylabel='Cumulative fraction')
    facet_title(ax, 'Rank distribution')

    # The MRR travels with the graph name in one legend under both panels.
    if len(kgs) > 1:
        fig.legend(*axes[1].get_legend_handles_labels(), loc='outside lower center',
                   ncol=3 if len(kgs) > 3 else len(kgs))
    panel_letters(fig, axes)
    return fig
