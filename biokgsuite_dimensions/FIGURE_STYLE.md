# Figure style

Every figure in the repository (dimension notebooks, the grounding experiment,
the supplementary embeddings and the set-construction tools) is drawn by a
function under `src/biokgsuite/viz/` on top of `viz/style.py`. The aim is a
Nature or Briefings in Bioinformatics display item: typeset at final size,
restrained, nothing decorative. Read `style.py` first; it is short.

## Rules

1. **No titles inside the figure.** The caption carries the title. The only
   heading a panel may carry is `facet_title(ax, text)` (7 pt, regular weight,
   left aligned), used when panels are facets of one thing (one per graph,
   per model, per metric).
2. **Panel letters** with `panel_label(ax, 'a')`: bold lowercase 8 pt, outside
   the top-left corner. Every multi-panel figure has them; a single panel has
   none.
3. **Typography** comes from `setup_style()`: 7 pt labels, 6 pt ticks and
   legends, Helvetica or Arial (Liberation Sans where those are absent). No
   fontsize above 8, nothing bold except panel letters, sentence case in every
   label ("Edge dropout (%)"), no raw keys in a label (`drug_gene` is
   "Drug–gene", `primekg` is "PrimeKG" via `KG_NAMES`).
4. **Colour follows the entity.** A graph is always `KG_PALETTE[kg]`
   (Okabe-Ito, ordered by `KG_ORDER`, distinguishable under the three common
   colour-vision deficiencies and in greyscale). Baselines, nulls, random pairs,
   the no-KG arm and chance lines are `NEUTRAL` or `FAINT` grey, dashed for
   reference lines. Heuristic scorers use `HEURISTIC_COLORS` (greys) because
   the scorer is not the subject. Magnitude uses `HEATMAP_CMAP`, one blue from
   light to dark; polarity uses `DIVERGING_CMAP` around a neutral midpoint.
   Text is never coloured by series.
5. **Axes**: `clean_ax(ax, ...)` on every axes: no gridlines, left and bottom
   hairline spines, outward ticks, horizontal tick labels (widen the figure or
   shorten the label rather than rotate). Large counts use `si_formatter()`
   ("1.5 M"), never an offset text. Percent axes say "(%)" with plain numbers.
6. **Marks**: bars 0.6 of the slot with no edge lines; a value printed on a bar
   only when the number is the message and there are at most about six bars;
   lines 1 pt with 3.5 pt markers; error bars 0.75 pt without caps; confidence
   bands at alpha 0.15 in the series colour; thresholds grey, dashed, with a
   6 pt muted label at the right end rather than a legend entry.
7. **Legends** frameless, 6 pt, placed clear of the data; direct labels instead
   of a legend when there are at most four series and room at the line end;
   never a legend for a single series.
8. **Size** from `fig_size(cols, rows)`: 89 mm single column or 183 mm double,
   about 50 mm per row, the smallest that keeps 6 pt text legible. Nothing
   overlaps; check the render.
9. **Six graphs.** Every function lays out for `len(KG_ORDER)` panels or bars
   and degrades to fewer. Facet grids wrap at three columns.
10. **Output** via `save_fig` (or `Outputs.figure`): PDF with editable text
    (fonttype 42) and a 600 dpi PNG, both from the same Figure.

## Checking a figure

Redraw from the checkpoint (`python tools/regenerate_figs.py --only NN`) or
from a small synthetic input, open the PNG, and confirm: nothing overlaps,
every label is a real name, panel letters are present, colours follow
`KG_PALETTE`, no text is larger than 8 pt, the width is 89 or 183 mm.
`tests/test_grounding_viz.py` renders the grounding figures on mock runs.
