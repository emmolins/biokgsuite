"""House style for every figure in the repository.

Typeset for print at final size, the way Nature-family and Briefings in
Bioinformatics artwork is prepared: Helvetica or Arial, 7 pt text with 6 pt
ticks (the journals' floor is 5 pt), 8 pt bold lowercase panel letters, 0.5 pt
hairline axes, no gridlines, no titles inside the figure (the caption carries
them), single-column 89 mm or double-column 183 mm widths, vector PDF with
embedded TrueType text plus a 600 dpi PNG.

Colour follows the entity. Every graph has one colour in every figure
(``KG_PALETTE``), taken from the Okabe-Ito set, which is distinguishable under
deuteranopia, protanopia and tritanopia and in greyscale print. Baselines,
nulls and reference lines are grey so they never compete with a graph.
Magnitude is one hue from light to dark (``HEATMAP_CMAP``); polarity is a
warm-cool pair around a neutral midpoint (``DIVERGING_CMAP``).
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ── Knowledge graphs: fixed order, one colour each (Okabe-Ito) ────────────────
KG_ORDER = ['primekg', 'hetionet', 'drkg', 'openbilink', 'biokg', 'matrix']
KG_NAMES = {'primekg': 'PrimeKG', 'hetionet': 'Hetionet', 'drkg': 'DRKG',
            'openbilink': 'OpenBioLink', 'biokg': 'BioKG', 'matrix': 'MATRIX'}
KG_PALETTE = {
    'primekg':    '#0072B2',   # blue
    'hetionet':   '#D55E00',   # vermilion
    'drkg':       '#009E73',   # green
    'openbilink': '#E69F00',   # orange
    'biokg':      '#CC79A7',   # purple
    'matrix':     '#56B4E9',   # sky blue
}
PALETTE = [KG_PALETTE[k] for k in KG_ORDER]

# ── Ink and reference colours ────────────────────────────────────────────────
TEXT_COLOR = '#1A1A1A'   # labels, values, legends: text never wears a series colour
TICK_COLOR = '#1A1A1A'
MUTED = '#6E6E6E'        # secondary text, annotations
FAINT = '#BDBDBD'        # reference lines, de-emphasised marks
GRID_COLOR = '#E6E6E6'
NEUTRAL = '#8C8C8C'      # baselines and nulls (no-KG arm, random pairs, degree-matched null)
ALERT_RED = '#D55E00'    # thresholds that mean "below this is a problem"

# ── Magnitude and polarity ───────────────────────────────────────────────────
HEATMAP_STOPS = ['#F7FBFF', '#DEEBF7', '#C6DBEF', '#9ECAE1', '#6BAED6', '#3182BD', '#08519C', '#08306B']
HEATMAP_CMAP = mcolors.LinearSegmentedColormap.from_list('single_hue_blue', HEATMAP_STOPS)
DIVERGING_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'blue_grey_vermilion', ['#0072B2', '#F0F0F0', '#D55E00'])

# Link-prediction scorers keep their own fixed hues (never reused for a graph).
HEURISTIC_COLORS = {'Adamic-Adar': '#1A1A1A', 'Common Neighbors': '#6E6E6E', 'Jaccard': '#BDBDBD'}

# Aliases kept so older call sites keep importing; new code uses the names above.
KG_OCEAN_BLUE, KG_RUST, KG_SAGE, KG_AMBER, KG_PLUM = PALETTE[:5]
OI_BLACK, OI_ORANGE, OI_SKY_BLUE, OI_GREEN = '#000000', '#E69F00', '#56B4E9', '#009E73'
OI_BLUE, OI_VERMILION, OI_PINK = '#0072B2', '#D55E00', '#CC79A7'

# ── Sizes (Nature column widths) ─────────────────────────────────────────────
_MM = 1 / 25.4
SINGLE_COL_W = 89 * _MM       # 3.50 in
DOUBLE_COL_W = 183 * _MM      # 7.20 in
ROW_H_STD = 50 * _MM          # 1.97 in per row


def fig_size(cols: int = 2, rows: int = 1, row_h: float = None) -> tuple:
    """(width, height) in inches: single column for ``cols=1``, double for 2."""
    return (SINGLE_COL_W if cols == 1 else DOUBLE_COL_W), (row_h or ROW_H_STD) * rows


def setup_style():
    """Apply the house style. Call once at the top of a notebook or script."""
    mpl.rcParams.update({
        'font.family':         'sans-serif',
        # Liberation Sans is metrically identical to Arial and present where
        # Helvetica and Arial are not; DejaVu is the last resort only.
        'font.sans-serif':     ['Helvetica Neue', 'Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
        'font.size':           7,
        'axes.titlesize':      7,
        'axes.titleweight':    'normal',
        'axes.labelsize':      7,
        'xtick.labelsize':     6,
        'ytick.labelsize':     6,
        'legend.fontsize':     6,
        'legend.title_fontsize': 6,
        'legend.frameon':      False,
        'legend.handlelength': 1.4,
        'legend.handletextpad': 0.5,
        'legend.columnspacing': 1.0,
        'legend.borderaxespad': 0.3,
        'figure.dpi':          150,
        'savefig.dpi':         600,
        'savefig.bbox':        'tight',
        'savefig.pad_inches':  0.02,
        'savefig.facecolor':   'white',
        'figure.facecolor':    'white',
        'axes.facecolor':      'white',
        'axes.edgecolor':      TEXT_COLOR,
        'axes.linewidth':      0.5,
        'axes.labelcolor':     TEXT_COLOR,
        'axes.labelpad':       3,
        'axes.titlepad':       4,
        'axes.spines.top':     False,
        'axes.spines.right':   False,
        'axes.grid':           False,
        'axes.axisbelow':      True,
        'axes.prop_cycle':     plt.cycler(color=PALETTE),
        'axes.formatter.use_mathtext': True,
        'xtick.color':         TEXT_COLOR,
        'ytick.color':         TEXT_COLOR,
        'xtick.direction':     'out',
        'ytick.direction':     'out',
        'xtick.major.size':    2.5,
        'ytick.major.size':    2.5,
        'xtick.major.width':   0.5,
        'ytick.major.width':   0.5,
        'xtick.major.pad':     2,
        'ytick.major.pad':     2,
        'lines.linewidth':     1.0,
        'lines.markersize':    3.5,
        'lines.markeredgewidth': 0.5,
        'patch.linewidth':     0,
        'errorbar.capsize':    0,
        'hatch.linewidth':     0.4,
        # Text stays text in the PDF, so a journal's art department can edit it.
        'pdf.fonttype':        42,
        'ps.fonttype':         42,
        'svg.fonttype':        'none',
    })


def setup_manuscript_style():
    """Same style; kept for scripts that called the older two-tier API."""
    setup_style()


def clean_ax(ax, title='', xlabel='', ylabel='', grid_axis='none'):
    """Open axes: left and bottom hairlines, outward ticks, no grid.

    ``title`` is drawn as a small facet label (7 pt, regular weight, left
    aligned), which is the only heading a panel carries; figure titles belong
    in the caption.
    """
    if title:
        facet_title(ax, title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(False)
    if grid_axis in ('x', 'y', 'both'):
        ax.grid(axis=grid_axis, color=GRID_COLOR, linewidth=0.4)
    ax.set_axisbelow(True)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(0.5)
        ax.spines[spine].set_color(TEXT_COLOR)
    ax.tick_params(axis='both', which='major', length=2.5, width=0.5, colors=TEXT_COLOR, labelsize=6)


def facet_title(ax, text, loc='left', pad=3):
    """A small, regular-weight heading naming what a panel shows (a graph, a model)."""
    ax.set_title(text, fontsize=7, fontweight='normal', color=TEXT_COLOR, loc=loc, pad=pad)


def panel_label(ax, label, x=-0.18, y=1.06, fontsize=8):
    """Bold lowercase panel letter outside the top-left corner of the axes."""
    ax.text(x, y, str(label).lower(), transform=ax.transAxes, fontsize=fontsize,
            fontweight='bold', color=TEXT_COLOR, va='bottom', ha='left')


def bar_labels(ax, bars, fmt='{:.2f}', offset=0.01, fontsize=6):
    """Value above each bar; use only when the number itself is the message."""
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, h + offset, fmt.format(h),
                ha='center', va='bottom', fontsize=fontsize, color=TEXT_COLOR)


def si_formatter(decimals=1):
    """Axis formatter writing 1500000 as 1.5 M and 2400 as 2.4 k, never 1e6."""
    def fmt(x, _pos):
        for div, suffix in ((1e9, ' G'), (1e6, ' M'), (1e3, ' k')):
            if abs(x) >= div:
                v = x / div
                return f'{v:.{decimals}f}{suffix}'.replace(f'.{"0" * decimals} ', ' ')
        return f'{x:g}'
    return FuncFormatter(fmt)


def kg_colour(kg: str) -> str:
    return KG_PALETTE.get(kg, NEUTRAL)


def kg_name(kg: str) -> str:
    return KG_NAMES.get(kg, kg)


def save_fig(fig, figs_dir, name, dpi=600):
    """Write ``name.pdf`` (vector, editable text) and ``name.png`` at ``dpi``."""
    figs_dir = Path(figs_dir)
    figs_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(figs_dir / f'{name}.{ext}', dpi=dpi, bbox_inches='tight', facecolor='white')
    print(f'  saved {name}.pdf and .png ({dpi} dpi)')
