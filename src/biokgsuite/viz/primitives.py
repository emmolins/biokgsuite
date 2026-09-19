"""Shared figure primitives for the dimension figures.

Everything here builds on ``biokgsuite.viz.style`` and adds only what more
than one dimension module needs: graph ordering and lookup, human-readable
labels for internal keys, reference lines with an end label, panel letters
that clear long tick labels, and the coverage heatmap that notebooks 01 and 02
share. Split out of the former eval_figs.py so each dimension owns only its own
figure code.
"""
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import (DOUBLE_COL_W, FAINT, HEATMAP_CMAP, KG_NAMES, KG_ORDER, MUTED,
                    ROW_H_STD, SINGLE_COL_W, TEXT_COLOR, facet_title, kg_colour,
                    panel_label)

# The fixed order and the published names live in style.py; these aliases keep
# the older import sites (dimension06, dimension07) working unchanged.
KGS = list(KG_ORDER)
NICE = dict(KG_NAMES)

# Entity words that appear in relation keys. A key made of two of them is a
# pair type ("drug_gene" is the drug-gene relation) and reads with an en dash;
# any other snake_case key is an attribute and reads as a sentence-case phrase.
_ENTITY_WORDS = {'drug', 'gene', 'protein', 'disease', 'pathway', 'phenotype',
                 'target', 'compound', 'symptom', 'anatomy', 'variant'}


def kgc(kg):
    """Colour of a graph; grey for anything that is not one of the six."""
    return kg_colour(kg)


def _present(cols):
    """The graphs in ``cols``, in ``KG_ORDER``; unknown keys are dropped."""
    return [k for k in KGS if k in list(cols)]


def human_label(key):
    """'drug_gene' -> 'Drug–gene', 'approval_year' -> 'Approval year'.

    Labels that are already prose ('Gene/Protein', 'Spearman r') pass through
    untouched: only snake_case keys are internal names that must not reach the
    manuscript.
    """
    s = str(key)
    if '_' not in s:
        return s
    parts = [p for p in s.split('_') if p]
    if parts and all(p.lower() in _ENTITY_WORDS for p in parts):
        return '–'.join(parts).capitalize()
    return ' '.join(parts).capitalize()


def truncate(text, limit=28):
    """Shorten a label to ``limit`` characters with an ellipsis."""
    s = str(text)
    return s if len(s) <= limit else s[:limit - 1].rstrip() + '…'


def _xlabels(ax, names):
    """Graph names on the x axis, horizontal; the figure is sized so they fit."""
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([NICE.get(n, n) for n in names])


def _title(ax, title, subtitle=None, pad=3):
    """Kept for older call sites: draws a facet title only. Figure titles and
    subtitles belong in the caption, so ``subtitle`` is accepted and ignored."""
    facet_title(ax, title, pad=pad)


def _df(data, index='KG'):
    df = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
    if index in df.columns:
        df = df.set_index(index)
    return df


def _save(fig, folder, names):
    os.makedirs(folder, exist_ok=True)
    names = [names] if isinstance(names, str) else names
    for nm in names:
        for e in ('png', 'pdf'):
            fig.savefig(os.path.join(folder, f'{nm}.{e}'), dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def facet_grid(n, ncol_max=3, row_h=None):
    """Figure and flat axes for ``n`` facets, wrapping at ``ncol_max`` columns.

    One facet fits a single column; anything wider takes the double column so
    that six graphs sit in two rows of three at the same panel size the paper
    prints. Unused axes are hidden, not deleted, so callers can still index by
    facet position.
    """
    ncol = min(ncol_max, max(1, n))
    nrow = int(np.ceil(n / ncol))
    width = SINGLE_COL_W if ncol == 1 else DOUBLE_COL_W
    fig, axes = plt.subplots(nrow, ncol, figsize=(width, (row_h or ROW_H_STD) * nrow),
                             squeeze=False, constrained_layout=True)
    axes = axes.reshape(-1)
    for ax in axes[n:]:
        ax.set_visible(False)
    return fig, axes, ncol


def ref_line(ax, y, label, color=FAINT, pad_pt=2):
    """Grey dashed reference line with a 6 pt muted label just past its right end.

    The label sits outside the axes, so it never lands on a bar or a marker
    that happens to reach the threshold, and it replaces a legend entry.
    """
    ax.axhline(y, color=color, ls='--', lw=0.6, zorder=1)
    ax.annotate(label, xy=(1, y), xycoords=('axes fraction', 'data'),
                xytext=(pad_pt, 0), textcoords='offset points',
                ha='left', va='center', fontsize=6, color=MUTED,
                annotation_clip=False)


def panel_letters(fig, axes, letters='abcdefghijklmnop', pad_pt=3):
    """Panel letters that clear each panel's tick labels.

    ``panel_label`` puts the letter a fixed fraction left of the axes, which
    lands on top of long y tick labels (source names, relation names). Drawing
    once gives the real extent of each panel's labels, so the letter can sit at
    the panel's left edge at the same height as the facet title.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax, letter in zip(axes, letters):
        if not ax.get_visible():
            continue
        full = ax.get_tightbbox(renderer)
        frame = ax.get_window_extent(renderer)
        px_per_pt = fig.dpi / 72
        # A panel without y tick labels (shared axis) has no room of its own,
        # so the letter still moves one letter width clear of the facet title.
        dx_px = min(full.x0 - frame.x0, -9 * px_per_pt)
        panel_label(ax, letter, x=dx_px / frame.width,
                    y=1 + pad_pt * px_per_pt / frame.height)


def _dark_cell(cmap, frac):
    """True when text on this cell must be white to stay legible."""
    r, g, b, _ = cmap(float(np.clip(frac, 0, 1)))
    return 0.299 * r + 0.587 * g + 0.114 * b < 0.5


# ───────────────────────── 01 / 02 · coverage heatmap ─────────────────────────
def coverage_heatmap(pivot, title, subtitle=None, not_scored=(), figsize=None,
                     cbar_label='Coverage (%)'):
    """Heatmap of percentages, rows x graphs, with explicit n/a cells.

    Entity and relation coverage (and entity validity) are the same picture
    over different rows, so they share one renderer: a difference between the
    figures would read as a difference in the measurement.

    NaN is drawn as a flat grey 'n/a' rather than as zero, because "not
    evaluable for this graph" and "evaluable and absent" are different findings.
    Rows named in ``not_scored`` are labelled as such, so a reader can see at a
    glance which rows feed the dimension score and which are reported beside it.

    ``title`` and ``subtitle`` are accepted for compatibility with the notebooks
    but not drawn: the caption carries them. ``cbar_label`` is likewise kept;
    the colourbar is labelled with the unit only, the caption says what it is.
    """
    # Graphs in the fixed order of the paper, whatever order the pivot arrived
    # in (a pandas pivot sorts its columns alphabetically).
    cols = _present(pivot.columns) + [c for c in pivot.columns if c not in KGS]
    pivot = pivot[cols]
    rows = list(pivot.index)
    data = pivot.values.astype(float)
    valid = ~np.isnan(data)

    # Up to four graph names fit horizontally in a single column; six need the
    # double column. Row height is fixed so a 3-row and a 6-row panel keep the
    # same cell size instead of one being squashed to a fixed canvas.
    if figsize is None:
        width = SINGLE_COL_W if len(cols) <= 4 else DOUBLE_COL_W
        figsize = (width, 0.28 * len(rows) + 0.3)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    im = ax.pcolormesh(np.where(valid, data, np.nan), cmap=HEATMAP_CMAP,
                       vmin=0, vmax=100, edgecolors='none', linewidth=0,
                       antialiased=False)
    ax.pcolormesh(np.where(~valid, 1, np.nan),
                  cmap=mcolors.ListedColormap(['#F2F2F2']), vmin=0, vmax=1,
                  edgecolors='none', linewidth=0, antialiased=False)

    # Published names, not the internal config keys: the key for OpenBioLink is
    # spelt 'openbilink' and would otherwise reach a figure in the manuscript.
    ax.set_xticks([x + 0.5 for x in range(len(cols))])
    ax.set_xticklabels([NICE.get(c, str(c)) for c in cols], fontsize=6, color=TEXT_COLOR)
    ax.set_yticks([y + 0.5 for y in range(len(rows))])
    ax.set_yticklabels([human_label(r) + ('\n(not scored)' if r in not_scored else '')
                        for r in rows], fontsize=6, color=TEXT_COLOR)
    ax.tick_params(length=0, which='both', pad=3)
    ax.grid(False, which='both')
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i in range(len(rows)):
        for j in range(len(cols)):
            val = data[i, j]
            if valid[i, j]:
                ax.text(j + 0.5, i + 0.5, f'{val:.1f}', ha='center', va='center',
                        fontsize=6,
                        color='white' if _dark_cell(HEATMAP_CMAP, val / 100) else TEXT_COLOR)
            else:
                ax.text(j + 0.5, i + 0.5, 'n/a', ha='center', va='center',
                        fontsize=6, color=MUTED, style='italic')

    # The bar lives in an inset so it is exactly as tall as the cells.
    cax = ax.inset_axes([1.03, 0, 0.035, 1])
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=6, length=2, width=0.5, colors=TEXT_COLOR)
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(TEXT_COLOR)
    cbar.set_label('(%)', fontsize=6, color=TEXT_COLOR, labelpad=3)
    return fig
