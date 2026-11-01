# 00 · Benchmark summary

**Computes.** Nothing of its own. It reads the seven dimension checkpoints and
assembles the sub-metric matrix (19 sub-metrics x the graphs in
`config/kgs.yaml`), the dimension and category aggregates, the three headline
figures, and three readings of whether the dimensions are redundant: Kendall's
*W* under five specifications, rank dispersion within each graph, and
dissociation counts. Run last.

**Inputs.** `config/kgs.yaml`, `results/checkpoints/01_coverage.pkl` through
`results/checkpoints/07_generalization.pkl`
**Outputs.** `results/{figures,tables}/00_summary/`
**Display item.** Figs 2 and 3; Supplementary Tables S3a and S4

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | nothing: every value is read from the checkpoint its dimension notebook wrote | the permutation *p*-value for Kendall's *W*, the 0.5 and 0.8 reference lines on the category chart, the red mark on summary-table values below 0.5 |

The only arithmetic performed here is the unweighted mean:

$$\text{Category}_c = \frac{1}{|D_c|}\sum_{d \in D_c} \text{score}_d, \qquad
\text{Overall} = \frac{1}{|D|}\sum_{d \in D} \text{score}_d$$

where $D$ is the set of dimensions the graph was actually scored on. A dimension
a graph was not scored on is absent from the sum rather than entering it as zero,
so the denominator is per graph and the table prints which dimensions each mean
rests on.

## Design notes

### The overall mean is over chance-floored dimensions

Every scored sub-metric runs from 0 at chance to 1. The AUROC-valued ones
(link prediction in 06, all three axes of 07) are mapped to 2·AUROC − 1, the
Gini coefficient, before their dimension is averaged; the proportion-valued
ones already start at 0. This puts the seven dimensions on one scale for the
radar, the heatmap and the overall mean, at the cost that the generalisation
dimension sits lower than a raw-AUROC reading would put it (a raw 0.78 is a
scored 0.55). Methods should state the map in one sentence, and the text
should quote raw AUROCs, which every 06 and 07 table keeps.

### A partial benchmark is refused rather than rendered

A missing checkpoint drops a spoke from the radar, a block from the heatmap and a
judge from the concordance, and the overall mean is then taken over a different
set of dimensions than the manuscript describes. All of it still renders, and
none of it says so on the figure. The loading cell therefore raises unless
`ALLOW_PARTIAL_DIMENSIONS` is set in the config cell.

The same argument applies one level down. Each checkpoint carries the graphs its
own notebook managed to score, so a dimension that saw two of six contributes two
rows here while the table around it looks complete. `assert_full_panel` checks
the union of scored graphs against `config/kgs.yaml` and reports, per dimension,
how many graphs it covered.

### N/A is not zero

An unmeasured quantity is carried as `np.nan` everywhere: the radar breaks the
polygon at a NaN vertex rather than closing it at the origin, the heatmap shades
the cell instead of colouring it on the score scale, and the means above skip it.
Zero is a score a graph can earn; not having been measured is not.

### Concordance is an effect size, not a *p*-value

Kendall's *W* is reported with its exact conversion to the mean pairwise Spearman
correlation, $\bar\rho = (kW - 1)/(k - 1)$. The permutation null permutes graph
labels within each metric, which preserves each metric's rank structure and
destroys any consistent ordering of graphs, but the metrics are not independent
judges: the three task-performance metrics share a scoring heuristic and the
coverage metrics share gold standards. The *p*-value is therefore optimistic, and
what the manuscript reports is the estimate and its stability across the five
specifications, not the *p*-value. Note also that $E[W] = 1/k$ under this null,
so *W* has a positive floor that rises as judges are removed, which is why the
jackknife and one-metric-per-dimension specifications carry their own null
column.

### What "discriminating" means

A metric enters the discriminating set when its range across the graphs exceeds
`_DISC_TOL = 0.05`, a floor on what counts as a real difference rather than
rounding, set a priori at 5% of the 0-1 scale every sub-metric is reported on.
It is a selection rule and not a threshold anything is scored against; the
jackknife and the exhaustive one-per-dimension enumeration show what the estimate
does when the set it selects is perturbed.

### Dissociation needs no null

At this panel size the assumption-free reading is the most defensible one. The
dissociation count asks only whether two properties order a pair of graphs the
same way, so it needs no distribution and no null: 0% means the two properties
are redundant, and 50% is what unrelated orderings give.

### No graph is emphasised

`HIGHLIGHT` is `None` on the radar. Drawing one series heavier or more opaque
than the rest puts a thumb on the comparison the benchmark exists to make, and
the graph being emphasised was the authors' own.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/00_summary/notebook.ipynb
```

Run after 01 through 07. The checkpoint reader and the results paths live in
`src/biokgsuite/io/outputs.py`, the panel guard in
`src/biokgsuite/io/config.py`, and the palette, colour map and figure sizes in
`src/biokgsuite/viz/style.py`.
