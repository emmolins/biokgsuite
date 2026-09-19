"""Figures for notebook 00, the benchmark summary (the manuscript's Figure 1).

The summary notebook measures nothing: every number here is read from the
checkpoint a dimension notebook wrote. The module therefore also owns the
three definitions that say what a summary is made of (``DIM_LABELS``,
``CATEGORIES`` and ``SUB_METRICS``) and the extraction helpers that turn a set
of checkpoints into figure inputs, so the notebook and
``tools/regenerate_figs.py`` read the checkpoints through one piece of code
and cannot drift apart.

Each function returns its Figure; saving goes through ``biokgsuite.io.Outputs``.
"""
import math

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

from biokgsuite.viz.primitives import NICE, _dark_cell, _present, kgc, ref_line
from biokgsuite.viz.style import (DOUBLE_COL_W, FAINT, GRID_COLOR, HEATMAP_CMAP, MUTED,
                                  SINGLE_COL_W, TEXT_COLOR, clean_ax, fig_size)

# ───────────────────────── what a summary is made of ────────────────────────────
DIM_LABELS = {
    'coverage':            'Coverage',
    'annotation_accuracy': 'Annotation accuracy',
    'trustworthiness':     'Trustworthiness',
    'topology':            'Topology',
    'stability':           'Stability',
    'task_performance':    'Task performance',
    'generalisation':      'Generalisation',
}

# Two-line forms for the labels that do not fit a radar spoke or a two-cell
# heatmap group at 6 pt; every other label is used as written.
TWO_LINE_LABELS = {
    'annotation_accuracy': 'Annotation\naccuracy',
    'task_performance':    'Task\nperformance',
}

CATEGORIES = {
    'Content':   ['coverage', 'annotation_accuracy', 'trustworthiness'],
    'Structure': ['topology', 'stability'],
    'Inference': ['task_performance', 'generalisation'],
}


def _sub(key):
    """Extractor reading one key of a checkpoint's ``sub_scores[kg]``."""
    return lambda ck, kg: ck.get('sub_scores', {}).get(kg, {}).get(key)


# Every sub-metric the seven dimension notebooks score, in their order, under
# the names their own title cells give them: (column label, dimension key,
# extractor(ck, kg)). A key renamed in a dimension notebook shows up as a
# blank column rather than as a silently different number.
SUB_METRICS = [
    ('Entity coverage',             'coverage',            _sub('entity_coverage')),
    ('Relation coverage',           'coverage',            _sub('relation_coverage')),
    ('Entity validity',             'annotation_accuracy', _sub('entity_validity')),
    ('Relational consistency',      'annotation_accuracy', _sub('relational_consistency')),
    ('Edge traceability',           'trustworthiness',     _sub('edge_traceability')),
    ('Source diversity',            'trustworthiness',     _sub('source_diversity')),
    ('Uncertainty quantification',  'trustworthiness',     _sub('uq_coverage')),
    ('Graph cohesion',              'topology',            _sub('connectedness')),
    ('Clustering coefficient',      'topology',            _sub('small_world')),
    ('Graph reachability',          'topology',            _sub('reachability')),
    ('Community purity',            'topology',            _sub('community_purity')),
    ('Random-dropout stability',    'stability',           _sub('random_stability')),
    ('Periphery-dropout stability', 'stability',           _sub('periphery_stability')),
    ('Link prediction',             'task_performance',    _sub('link_prediction')),
    ('Neighbourhood retrieval',     'task_performance',    _sub('nbhd_retrieval')),
    ('Multi-hop reasoning',         'task_performance',    _sub('multihop_reasoning')),
    ('Data-sparse generalisation',  'generalisation',      _sub('data_sparse_gen')),
    ('Cross-domain generalisation', 'generalisation',      _sub('cross_domain_gen')),
    ('Prospective generalisation',  'generalisation',      _sub('temporal_gen')),
]


# ───────────────────────── checkpoint extraction ────────────────────────────────
def _is_num(v):
    return v is not None and not (isinstance(v, float) and math.isnan(v))


def dimension_scores(ckpts, kgs):
    """``summary_scores[kg][dim]`` from ``{dim: checkpoint}``, on 0-1.

    A dimension a graph was not scored on is absent rather than 0.0, so it
    stays out of every mean instead of being averaged in as a failure.
    """
    scores = {kg: {} for kg in kgs}
    for dim, ck in ckpts.items():
        for kg in kgs:
            v = ck.get('summary_scores', {}).get(kg)
            if v is not None:
                scores[kg][dim] = float(v)
    return scores


def category_means(summary_scores, categories=None):
    """``category_scores[kg][cat]``: mean over the loaded dimensions of a category."""
    categories = categories or CATEGORIES
    out = {}
    for kg, dims in summary_scores.items():
        out[kg] = {}
        for cat, cat_dims in categories.items():
            vals = [dims[d] for d in cat_dims if d in dims]
            out[kg][cat] = round(float(np.mean(vals)), 4) if vals else None
    return out


def metric_matrix(ckpts, kgs, sub_metrics=None, dim_labels=None):
    """Inputs of the metric heatmap, read straight from the checkpoints.

    Returns ``(heat_data, sub_metric_labels, dim_of_metric, dim_data,
    dim_col_labels)``. ``heat_data`` is graphs x sub-metrics; ``dim_data`` is
    graphs x (dimensions + overall) and is recomputed as the unweighted mean of
    the sub-metric cells, so the block beside the heatmap is arithmetic the
    reader can check against the cells to its left.
    """
    sub_metrics = sub_metrics or SUB_METRICS
    dim_labels = dim_labels or DIM_LABELS
    heat = np.full((len(kgs), len(sub_metrics)), np.nan)
    for j, (_, dim, extract) in enumerate(sub_metrics):
        ck = ckpts.get(dim)
        if ck is None:
            continue
        for i, kg in enumerate(kgs):
            v = extract(ck, kg)
            if _is_num(v):
                heat[i, j] = float(v)

    dim_of_metric = [d for _, d, _ in sub_metrics]
    dim_order = list(dict.fromkeys(dim_of_metric))
    dim_data = np.full((len(kgs), len(dim_order) + 1), np.nan)
    for i in range(len(kgs)):
        for k, d in enumerate(dim_order):
            vals = heat[i, [j for j, dd in enumerate(dim_of_metric) if dd == d]]
            vals = vals[~np.isnan(vals)]
            if vals.size:
                dim_data[i, k] = vals.mean()
        dims = dim_data[i, :len(dim_order)]
        dims = dims[~np.isnan(dims)]
        if dims.size:
            dim_data[i, -1] = dims.mean()
    labels = [m[0] for m in sub_metrics]
    col_labels = [dim_labels.get(d, d) for d in dim_order] + ['Overall']
    return heat, labels, dim_of_metric, dim_data, col_labels


# ───────────────────────── shared drawing helpers ───────────────────────────────
def _kg_order(kgs):
    """The requested graphs in ``KG_ORDER``, then any unknown key as given."""
    kgs = list(kgs)
    return _present(kgs) + [k for k in kgs if k not in _present(kgs)]


def _legend_below(fig, kgs, ncol=3, marker=False, anchor=None):
    """6 pt frameless legend under the panel, one entry per graph.

    With ``anchor`` (figure fraction) the legend hangs from that point, for a
    figure laid out by hand; otherwise constrained layout makes room for it.
    """
    if marker:
        handles = [Line2D([0], [0], marker='o', ls='none', ms=3.5, color=kgc(k),
                          mec='white', mew=0.4, label=NICE.get(k, k)) for k in kgs]
    else:
        handles = [Line2D([0], [0], color=kgc(k), lw=1, label=NICE.get(k, k)) for k in kgs]
    if anchor is None:
        fig.legend(handles=handles, loc='outside lower center', ncol=min(ncol, len(kgs)),
                   fontsize=6, frameon=False)
    else:
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=anchor,
                   ncol=min(ncol, len(kgs)), fontsize=6, frameon=False)


def _overlaps(a, b, pad=1.0):
    """True when two display-space bounding boxes intersect (``pad`` px slack)."""
    return not (a.x1 + pad < b.x0 or b.x1 + pad < a.x0 or a.y1 + pad < b.y0 or b.y1 + pad < a.y0)


# ───────────────────────── 00 · dimension scores (Fig 1) ────────────────────────
def fig00_dimension_scores(summary_scores, dims, dim_labels, kgs):
    """Dot plot of every dimension score, one row per dimension, one dot per graph.

    The bottom row is the unweighted mean of the dimensions a graph was scored
    on; a missing dimension is skipped, never counted as zero. Graph names are
    written beside the dots of that row when they fit without touching another
    dot or name; otherwise a legend names them.
    """
    kgs = _kg_order(kgs)
    dims = list(dims)
    n = len(dims)
    gap = 0.7
    y_overall = n + gap
    row_y = list(range(n)) + [y_overall]

    fig, ax = plt.subplots(figsize=fig_size(1, 1, row_h=70 / 25.4), constrained_layout=True)

    # Per-row values: dimension scores, then the overall mean per graph.
    rows = []
    for d in dims:
        rows.append({kg: summary_scores.get(kg, {}).get(d, np.nan) for kg in kgs})
    overall = {}
    for kg in kgs:
        vals = [summary_scores[kg][d] for d in dims if _is_num(summary_scores.get(kg, {}).get(d))]
        overall[kg] = float(np.mean(vals)) if vals else np.nan
    rows.append(overall)

    for y, vals in zip(row_y, rows):
        finite = [v for v in vals.values() if _is_num(v)]
        if not finite:
            ax.text(0.5, y, 'n/a', ha='center', va='center', fontsize=6, color=MUTED,
                    style='italic', zorder=3)
            continue
        # A reference tick at the midpoint of the scale, so a row reads
        # "above or below one half" without a gridline.
        ax.plot([0.5, 0.5], [y - 0.28, y + 0.28], color=FAINT, lw=0.5, zorder=1,
                solid_capstyle='butt')
        if len(finite) > 1:
            ax.plot([min(finite), max(finite)], [y, y], color=FAINT, lw=0.6, zorder=2,
                    solid_capstyle='butt')
        for kg, v in vals.items():
            if _is_num(v):
                ax.plot(v, y, marker='o', ms=4, ls='none', color=kgc(kg), mec='white',
                        mew=0.4, zorder=4 if y == y_overall else 3)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(y_overall + 0.6, -0.6)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
    ax.set_yticks(row_y)
    ax.set_yticklabels([dim_labels.get(d, d) for d in dims] + ['Overall'])
    clean_ax(ax, xlabel='Dimension score')
    ax.tick_params(axis='y', length=0)
    ax.spines['left'].set_visible(False)

    # Direct labels beside the overall dots, or a legend when they collide.
    # Each label is anchored to its dot in data coordinates with a point
    # offset, so it survives the dpi change at save time; the collision test
    # uses the renderer, so it is made on the real text width.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_per_pt = fig.dpi / 72
    fig_box = fig.bbox
    dot_r = 2.2 * px_per_pt
    dots = {kg: v for kg, v in overall.items() if _is_num(v)}
    taken = []
    for v in dots.values():
        x, yy = ax.transData.transform((v, y_overall))
        taken.append(Bbox([[x - dot_r, yy - dot_r], [x + dot_r, yy + dot_r]]))
    if len(dots) > 1:        # the connector line between the outer dots
        x0, yy = ax.transData.transform((min(dots.values()), y_overall))
        x1, _ = ax.transData.transform((max(dots.values()), y_overall))
        taken.append(Bbox([[x0, yy - 1], [x1, yy + 1]]))
    taken.append(Bbox([ax.transData.transform((0.5, y_overall + 0.28)),
                       ax.transData.transform((0.5, y_overall - 0.28))]))
    placed = []
    ok = bool(dots)
    for kg, v in sorted(dots.items(), key=lambda kv: kv[1]):
        name = NICE.get(kg, kg)
        found = None
        # Beside the dot when the row has room, otherwise just above it (centred,
        # or ending or starting at the dot near the edges), in the gap that
        # separates the overall row from the dimensions.
        candidates = (((4.2, 0), 'left', 'center'), ((-4.2, 0), 'right', 'center'),
                      ((0, 4.2), 'center', 'bottom'), ((2, 4.2), 'right', 'bottom'),
                      ((-2, 4.2), 'left', 'bottom'))
        for xytext, ha, va in candidates:
            t = ax.annotate(name, xy=(v, y_overall), xytext=xytext,
                            textcoords='offset points', fontsize=6, color=TEXT_COLOR,
                            ha=ha, va=va, zorder=5, annotation_clip=False)
            bb = t.get_window_extent(renderer)
            # The label may run into the figure margin, but not past the
            # 89 mm edge, which is where the saved figure ends.
            inside = fig_box.x0 <= bb.x0 and bb.x1 <= fig_box.x1
            if inside and not any(_overlaps(bb, o) for o in taken):
                found = t
                taken.append(bb)
                break
            t.remove()
        if found is None:
            ok = False
            break
        placed.append(found)
    if not ok:
        for t in placed:
            t.remove()
        _legend_below(fig, kgs, ncol=3, marker=True)
    return fig


# ───────────────────────── 00 · quality profiles (radar) ────────────────────────
def fig00_quality_profiles(summary_scores, dims, dim_labels, kgs):
    """One polygon per graph over the dimension scores.

    A dimension a graph was not scored on is a NaN vertex: matplotlib breaks
    the polygon there, so an unmeasured dimension is never drawn at the origin
    as a worst-possible score. No graph is emphasised over the rest.
    """
    kgs = _kg_order(kgs)
    dims = list(dims)
    n = len(dims)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    closed = angles + angles[:1]

    # An 89 mm square. The polar axes is placed by hand (below) because
    # constrained layout centres a polar panel in its slot and leaves a band
    # of white between the lowest spoke labels and the legend.
    fig, ax = plt.subplots(figsize=(SINGLE_COL_W, SINGLE_COL_W),
                           subplot_kw={'polar': True})
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=6, color=MUTED)
    # Radial numbers sit in the gap just before the first spoke, with a white
    # halo so a polygon edge passing underneath never cuts through them.
    ax.set_rlabel_position(np.degrees(angles[0] - np.pi / n) if n else 0)
    for lab in ax.get_yticklabels():
        lab.set_bbox(dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.8))
    ax.spines['polar'].set_visible(False)
    ax.grid(color=GRID_COLOR, lw=0.4)

    for kg in kgs:
        vals = [summary_scores.get(kg, {}).get(d, np.nan) for d in dims]
        vals = [v if _is_num(v) else np.nan for v in vals]
        ring = vals + vals[:1]
        ax.plot(closed, ring, color=kgc(kg), lw=1, zorder=2, label=NICE.get(kg, kg))
        ax.fill(closed, ring, color=kgc(kg), alpha=0.06, zorder=1, lw=0)

    ax.set_xticks(angles)
    ax.set_xticklabels([TWO_LINE_LABELS.get(d, dim_labels.get(d, d)) for d in dims],
                       fontsize=6, color=TEXT_COLOR)
    ax.tick_params(axis='x', pad=2)
    for lab, ang in zip(ax.get_xticklabels(), angles):
        sx, sy = np.cos(np.pi / 2 - ang), np.sin(np.pi / 2 - ang)
        lab.set_horizontalalignment('left' if sx > 0.1 else 'right' if sx < -0.1 else 'center')
        lab.set_verticalalignment('bottom' if sy > 0.1 else 'top' if sy < -0.1 else 'center')

    # Largest circle whose spoke labels still fit the width: start from the
    # full width, draw, and shrink by whatever the labels overhang. The legend
    # then sits just under the lowest label, so no band of white opens up.
    fig_w, fig_h = fig.get_size_inches()
    renderer = fig.canvas.get_renderer()
    d = 0.9 * fig_w
    for _ in range(3):
        ax.set_position([(fig_w - d) / 2 / fig_w, 1 - 0.02 - d / fig_h, d / fig_w, d / fig_h])
        fig.canvas.draw()
        bb = ax.get_tightbbox(renderer).transformed(fig.transFigure.inverted())
        over = max(-bb.x0, bb.x1 - 1, -0.0)
        if over <= 0.002:
            break
        d -= 2 * over * fig_w + 0.02
    _legend_below(fig, kgs, ncol=3, anchor=(0.5, bb.y0 - 0.015))
    return fig


# ───────────────────────── 00 · category scores ─────────────────────────────────
def fig00_category_scores(category_scores, categories, kgs):
    """Grouped bars: one group per category, one bar per graph.

    The lines at 0.5 and 0.8 are reading aids, not thresholds any score is
    judged against; nothing downstream uses them.
    """
    kgs = _kg_order(kgs)
    cats = list(categories)
    x = np.arange(len(cats))
    n_kg = max(len(kgs), 1)
    w = 0.6 / n_kg
    n_bars = len(cats) * n_kg
    # Three categories by six graphs (18 bars) still fit a single column;
    # more than that takes the double.
    fig, ax = plt.subplots(figsize=fig_size(1 if n_bars <= 18 else 2, 1),
                           constrained_layout=True)
    for i, kg in enumerate(kgs):
        vals = [category_scores.get(kg, {}).get(c) for c in cats]
        vals = [v if _is_num(v) else np.nan for v in vals]
        ax.bar(x + (i - (n_kg - 1) / 2) * w, vals, w, color=kgc(kg), zorder=2,
               label=NICE.get(kg, kg))
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_xlim(-0.5, len(cats) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ref_line(ax, 0.5, '0.5')
    ref_line(ax, 0.8, '0.8')
    clean_ax(ax, ylabel='Mean dimension score')
    _legend_below(fig, kgs, ncol=3)
    return fig


# ───────────────────────── 00 · metric and dimension heatmap ────────────────────
def _cells(ax, x0, data, vmin=0, vmax=1):
    """One block of cells starting at column ``x0``; returns the mappable."""
    n_r, n_c = data.shape
    valid = ~np.isnan(data)
    xe = np.arange(n_c + 1) + x0
    ye = np.arange(n_r + 1)
    im = ax.pcolormesh(xe, ye, np.where(valid, data, np.nan), cmap=HEATMAP_CMAP,
                       vmin=vmin, vmax=vmax, edgecolors='none', linewidth=0, antialiased=False)
    ax.pcolormesh(xe, ye, np.where(~valid, 1, np.nan),
                  cmap=mcolors.ListedColormap(['#F2F2F2']), vmin=0, vmax=1,
                  edgecolors='none', linewidth=0, antialiased=False)
    for i in range(n_r):
        for j in range(n_c):
            v = data[i, j]
            if valid[i, j]:
                ax.text(x0 + j + 0.5, i + 0.5, f'{v:.2f}', ha='center', va='center',
                        fontsize=6,
                        color='white' if _dark_cell(HEATMAP_CMAP, (v - vmin) / (vmax - vmin))
                        else TEXT_COLOR)
            else:
                ax.text(x0 + j + 0.5, i + 0.5, 'n/a', ha='center', va='center',
                        fontsize=6, color=MUTED, style='italic')
    return im


def fig00_metric_heatmap(heat_data, sub_metric_labels, dim_of_metric, dim_labels,
                         dim_data, dim_col_labels, kgs):
    """Every sub-metric per graph, with the dimension means and their overall
    mean set off to the right.

    An n/a cell is flat grey rather than a colour on the score scale, because
    "not evaluable for this graph" and "scored zero" are different findings.
    Column labels are rotated here only because 17 sub-metric names cannot
    stand upright in a double column.
    """
    # Rows are labelled in the order the matrices arrived in: the caller built
    # them, so reordering the names here would mislabel the data.
    kgs = list(kgs)
    heat = np.asarray(heat_data, dtype=float)
    dim = np.asarray(dim_data, dtype=float)
    n_rows, n_metric = heat.shape
    n_dim = dim.shape[1]
    gap = 0.6
    dim_off = n_metric + gap
    n_total = dim_off + n_dim

    # The cell is square and sized by the width, so the height follows the
    # number of graphs plus the space the rotated labels and headers take.
    cell_in = (DOUBLE_COL_W - 1.05) / n_total
    fig_h = cell_in * n_rows + 1.55
    fig, ax = plt.subplots(figsize=(DOUBLE_COL_W, fig_h), constrained_layout=True)

    im = _cells(ax, 0, heat)
    _cells(ax, dim_off, dim)
    ax.set_xlim(0, n_total)
    ax.set_ylim(n_rows, 0)
    ax.set_aspect('equal')

    # Column labels, rotated (the one permitted exception), and graph names.
    xt = [j + 0.5 for j in range(n_metric)] + [dim_off + j + 0.5 for j in range(n_dim)]
    ax.set_xticks(xt)
    ax.set_xticklabels(list(sub_metric_labels) + list(dim_col_labels), rotation=45,
                       ha='right', rotation_mode='anchor', fontsize=6, color=TEXT_COLOR)
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    ax.set_yticklabels([NICE.get(k, k) for k in kgs], fontsize=6, color=TEXT_COLOR)
    ax.tick_params(length=0, pad=3)
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

    # Group headers above the metric block: a faint bracket spanning the
    # group's columns and its dimension name above, 6 pt regular.
    groups = []
    for j, d in enumerate(dim_of_metric):
        if groups and groups[-1][0] == d:
            groups[-1][2] = j + 1
        else:
            groups.append([d, j, j + 1])
    groups = [(d, s, e) for d, s, e in groups] + [('__dims__', dim_off, dim_off + n_dim)]
    y_line, y_text = -0.3, -0.5
    for d, s, e in groups:
        ax.plot([s + 0.12, e - 0.12], [y_line, y_line], color=FAINT, lw=0.5,
                clip_on=False, solid_capstyle='butt')
        label = 'Dimension scores' if d == '__dims__' else dim_labels.get(d, d)
        t = ax.text((s + e) / 2, y_text, label, ha='center', va='bottom', fontsize=6,
                    color=TEXT_COLOR, clip_on=False)
        # A name wider than its bracket wraps onto two lines rather than
        # running into the neighbouring group's name.
        fig.canvas.draw()
        span = ax.transData.transform((e, 0))[0] - ax.transData.transform((s, 0))[0]
        if t.get_window_extent(fig.canvas.get_renderer()).width > span and ' ' in label:
            t.set_text(TWO_LINE_LABELS.get(d, label.replace(' ', '\n', 1)))

    # Colourbar as tall as the cells, just past the dimension block.
    cax = ax.inset_axes([1 + 0.25 / n_total, 0, 0.35 / n_total, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=6, length=2, width=0.5, colors=TEXT_COLOR)
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(TEXT_COLOR)
    cbar.set_label('Score', fontsize=6, color=TEXT_COLOR, labelpad=3)
    return fig
