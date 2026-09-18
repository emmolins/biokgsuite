"""Publication-ready plotting utilities for BioKGSuite.

Palette: Custom colorblind-safe categorical + warm-grey heatmap palette.
Dimensions: single-column 89 mm, double-column 183 mm (Nature standard).
DPI: 300 for print; 150 for screen preview.
Font: minimum 8 pt; titles 10 pt bold.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from pathlib import Path

# ── Categorical palette — knowledge graphs ────────────────────────────────────
# Colorblind-safe (deuteranopia + protanopia), equal visual weight on white.
KG_OCEAN_BLUE = '#2E6B8A'   # KG 1
KG_RUST       = '#C46A52'   # KG 2
KG_SAGE       = '#6A8A6B'   # KG 3
KG_AMBER      = '#C49A3C'   # KG 4 — warm amber, distinct from existing palette
KG_PLUM       = '#8B5E8B'   # KG 5 — muted plum, distinct from all above

# Text / axis colours
TEXT_COLOR  = '#333333'   # dark grey for titles, labels, annotations
TICK_COLOR  = '#333333'   # same dark grey for tick labels
ALERT_RED   = '#C46A52'   # rust — used for thresholds/warnings
GRID_COLOR  = '#E0E0E0'   # light grey grid lines

# ── Ordered palette for multi-series plots ───────────────────────────────────
PALETTE = [KG_OCEAN_BLUE, KG_RUST, KG_SAGE, KG_AMBER, KG_PLUM]

# ── Per-KG colours (consistent across all notebooks) ─────────────────────────
KG_PALETTE = {
    'primekg':    KG_OCEAN_BLUE,   # ocean blue
    'hetionet':   KG_RUST,          # rust
    'drkg':       KG_SAGE,          # sage
    'openbilink': KG_AMBER,         # amber
    'biokg':      KG_PLUM,          # plum
}

# ── Sequential heatmap palette — warm grey (7 stops) ─────────────────────────
HEATMAP_STOPS = [
    '#EDEADF', '#CCC7BA', '#A8A296', '#868074',
    '#655F55', '#48443C', '#2E2B26',
]
HEATMAP_CMAP = mcolors.LinearSegmentedColormap.from_list('warm_grey', HEATMAP_STOPS)

# Heuristic colours for link-prediction figures
HEURISTIC_COLORS = {
    'Common Neighbors': KG_OCEAN_BLUE,
    'Jaccard':          KG_RUST,
    'Adamic-Adar':      KG_SAGE,
}

# ── Backward-compatible aliases from old Okabe-Ito palette ───────────────────
OI_BLACK     = '#000000'
OI_ORANGE    = KG_RUST
OI_SKY_BLUE  = KG_OCEAN_BLUE
OI_GREEN     = KG_SAGE
OI_BLUE      = KG_OCEAN_BLUE
OI_VERMILION = KG_RUST
OI_PINK      = KG_SAGE

# ── Figure dimensions (Nature column widths) ──────────────────────────────────
_MM = 1 / 25.4                    # mm → inches
SINGLE_COL_W = 89  * _MM          # ~3.50 in
DOUBLE_COL_W = 183 * _MM          # ~7.20 in
ROW_H_STD    = 55  * _MM          # ~2.17 in per row (default)


def fig_size(cols: int = 2, rows: int = 1, row_h: float = None) -> tuple:
    """Return (width_in, height_in) for a Nature-standard figure.

    Parameters
    ----------
    cols : int
        1 → single-column (89 mm); 2 → double-column (183 mm).
    rows : int
        Number of subplot rows; height = rows × row_h.
    row_h : float, optional
        Row height in inches.  Defaults to ROW_H_STD (~2.2 in).
    """
    w = SINGLE_COL_W if cols == 1 else DOUBLE_COL_W
    h = (row_h or ROW_H_STD) * rows
    return w, h


def setup_style():
    """Apply the shared matplotlib style.  Call once at the top of each notebook."""
    mpl.rcParams.update({
        # Font
        'font.size':           9,
        'font.family':         'sans-serif',
        # Liberation Sans is metrically identical to Arial and is present on
        # machines that lack Helvetica/Arial; without it matplotlib falls
        # back to DejaVu Sans, which is wider and not a Nature-preferred face.
        'font.sans-serif':     ['Helvetica Neue', 'Helvetica', 'Arial',
                                'Liberation Sans', 'DejaVu Sans'],
        'axes.titlesize':      10,
        'axes.titleweight':    'bold',
        'axes.labelsize':      9,
        'xtick.labelsize':     8,
        'ytick.labelsize':     8,
        'legend.fontsize':     8,
        'legend.framealpha':   0,
        'legend.edgecolor':    'none',
        'legend.frameon':      False,
        # Resolution
        'figure.dpi':          150,
        'savefig.dpi':         300,
        'savefig.bbox':        'tight',
        'savefig.facecolor':   'white',
        # Axes
        'axes.facecolor':      'white',
        'figure.facecolor':    'white',
        'axes.titlepad':       8,
        'axes.labelpad':       4,
        'axes.spines.top':     False,
        'axes.spines.right':   False,
        'axes.edgecolor':      '#333333',
        'axes.linewidth':      0.5,
        'axes.prop_cycle':     plt.cycler(color=PALETTE),
        # Grid
        'axes.grid':           True,
        'grid.color':          GRID_COLOR,
        'grid.linewidth':      0.5,
        'grid.alpha':          1.0,
        'axes.axisbelow':      True,
        # Lines / patches
        'lines.linewidth':     1.5,
        'patch.linewidth':     0.5,
    })


def clean_ax(ax, title='', xlabel='', ylabel='', grid_axis='y'):
    """Apply consistent axis styling: grid, minimal spines, correct font sizes.

    Parameters
    ----------
    grid_axis : {'x', 'y', 'both', 'none'}
        Which axis to show grid lines on.
    """
    if title:
        ax.set_title(title, fontsize=10, fontweight='bold', color=TEXT_COLOR, pad=8)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9, color=TICK_COLOR, labelpad=4)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9, color=TICK_COLOR, labelpad=4)
    ax.tick_params(axis='both', labelsize=8, colors=TICK_COLOR)
    ax.tick_params(axis='y', length=0)
    ax.tick_params(axis='x', length=3, color='#333333')
    # Grid
    if grid_axis != 'none':
        ax.grid(axis=grid_axis, color=GRID_COLOR, linewidth=0.5, alpha=1.0)
    ax.set_axisbelow(True)
    # Spines
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color('#333333')
    ax.spines['bottom'].set_linewidth(0.5)


def panel_label(ax, label, x=-0.12, y=1.05, fontsize=11):
    """Add a bold panel letter (a, b, c … or A, B, C …) to axes.

    Parameters
    ----------
    ax : matplotlib Axes
    label : str — e.g. 'a', 'b', 'A', 'B'
    """
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', color=TEXT_COLOR,
            va='top', ha='right')


def bar_labels(ax, bars, fmt='{:.3f}', offset=0.02, fontsize=8):
    """Add value labels above bars."""
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            h + offset,
            fmt.format(h),
            ha='center', va='bottom',
            fontsize=fontsize, color=TEXT_COLOR,
        )


def setup_manuscript_style():
    """Nature-submission variant of `setup_style`.

    Same palette, spines and grid as the notebook style, but typeset for
    print at final size: 7 pt body (Nature's floor is 5 pt), hairline axes,
    and TrueType embedding so the PDF stays editable in Illustrator.

    Call INSTEAD of `setup_style()` in manuscript/figures/*.py. Notebook
    figures should keep `setup_style()` — 7 pt is unreadable on screen.
    """
    setup_style()
    mpl.rcParams.update({
        'font.size':        7,
        'axes.titlesize':   7.5,
        'axes.labelsize':   7,
        'xtick.labelsize':  6.5,
        'ytick.labelsize':  6.5,
        'legend.fontsize':  6.5,
        'axes.linewidth':   0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size':  2.5,
        'ytick.major.size':  2.5,
        'lines.linewidth':   1.2,
        # keep vector text as text, not outlines
        'pdf.fonttype':     42,
        'ps.fonttype':      42,
        'svg.fonttype':     'none',
    })


def save_fig(fig, figs_dir, name, dpi=300):
    """Save a figure as both PDF and PNG for publication.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    figs_dir : pathlib.Path
        Directory to save into (created if needed).
    name : str
        Base filename without extension, e.g. '01_entity_coverage'.
    dpi : int
        Raster resolution for the PNG (the PDF is vector, so this only
        affects embedded images). 300 for notebooks, 600 for manuscript
        submission.
    """
    figs_dir = Path(figs_dir)
    figs_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(figs_dir / f'{name}.{ext}', dpi=dpi,
                    bbox_inches='tight', facecolor='white')
    print(f'  → Saved: {name}.pdf / .png  ({dpi} dpi)')

