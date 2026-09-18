"""Shared figure primitives for the dimension figures.

Colour lookup, KG ordering/labels, title styling and multi-format saving.
Split out of the former eval_figs.py so each dimension owns only its own
figure code.
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from .style import setup_style, KG_PALETTE, TEXT_COLOR, GRID_COLOR, HEATMAP_CMAP, ALERT_RED

KGS = ['primekg', 'hetionet', 'drkg', 'openbilink', 'biokg', 'matrix']
NICE = {'primekg': 'PrimeKG', 'hetionet': 'Hetionet', 'drkg': 'DRKG',
        'openbilink': 'OpenBioLink', 'biokg': 'BioKG', 'matrix': 'Matrix'}


def kgc(kg):
    return KG_PALETTE.get(kg, '#888888')          # Matrix (+any extra) -> grey, as in nb 00


def _present(cols):
    return [k for k in KGS if k in list(cols)]


def _xlabels(ax, names):
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([NICE.get(n, n) for n in names], rotation=22, ha='right')


def _title(ax, title, subtitle=None, pad=10):
    ax.set_title(title, fontsize=11, fontweight='bold', loc='left', color=TEXT_COLOR, pad=pad)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=8.5, color='#777', va='bottom', ha='left')


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


# ───────────────────────── 02 · entity validity heatmap ─────────────────────────


# ───────────────────────── 01 · coverage heatmaps ─────────────────────────
def coverage_heatmap(pivot, title, subtitle=None, not_scored=(), figsize=None,
                     cbar_label='Coverage (%)'):
    """Heatmap of coverage percentages, rows x KGs, with explicit N/A cells.

    Entity and relation coverage are the same picture over different rows, so
    they share one renderer: a difference between the two figures would read as
    a difference in the measurement.

    NaN is drawn as a flat grey 'N/A' rather than as zero, because "not
    evaluable for this graph" and "evaluable and absent" are different findings.
    Rows named in ``not_scored`` are labelled as such, so a reader can see at a
    glance which rows feed the dimension score and which are reported beside it.
    """
    import matplotlib.colors as mcolors

    # Height follows the row count so a 3-row and a 6-row panel keep the same
    # cell size, rather than one being squashed to fit a fixed canvas.
    if figsize is None:
        figsize = (7.4, len(pivot.index) * 0.62 + 1.8)
    fig, ax = plt.subplots(figsize=figsize)
    data = pivot.values.astype(float)
    valid = ~np.isnan(data)

    im = ax.pcolormesh(np.where(valid, data, np.nan), cmap=HEATMAP_CMAP,
                       vmin=0, vmax=100, edgecolors='none', linewidth=0,
                       antialiased=False)
    ax.pcolormesh(np.where(~valid, 1, np.nan),
                  cmap=mcolors.ListedColormap(['#f0f0f0']), vmin=0, vmax=1,
                  edgecolors='none', linewidth=0, antialiased=False)

    ax.set_xticks([x + 0.5 for x in range(len(pivot.columns))])
    # Published names, not the internal config keys: the key for OpenBioLink is
    # spelt 'openbilink' and would otherwise reach a figure in the manuscript.
    ax.set_xticklabels([NICE.get(c, str(c)) for c in pivot.columns],
                       fontsize=9, color=TEXT_COLOR, rotation=22, ha='right')
    ax.set_yticks([y + 0.5 for y in range(len(pivot.index))])
    ax.set_yticklabels([f'{r}\n(not scored)' if r in not_scored else str(r)
                        for r in pivot.index], fontsize=9, color=TEXT_COLOR)
    ax.tick_params(length=0, which='both')
    ax.grid(False, which='both')
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i, row in enumerate(pivot.index):
        for j, col in enumerate(pivot.columns):
            val = pivot.loc[row, col]
            if pd.notna(val):
                ax.text(j + 0.5, i + 0.5, f'{val:.1f}%', ha='center', va='center',
                        fontsize=9, fontweight='bold',
                        color='white' if val >= 55 else TEXT_COLOR)
            else:
                ax.text(j + 0.5, i + 0.5, 'N/A', ha='center', va='center',
                        fontsize=8.5, color='#888888', style='italic')

    cbar = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(cbar_label, fontsize=8.5, color=TEXT_COLOR)
    _title(ax, title, subtitle, pad=34 if subtitle else 12)
    plt.tight_layout()
    return fig
