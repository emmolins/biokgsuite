# G03 · Analysis

**Computes.** For each model, the MRR gain of every graph over the `none` arm
with a paired bootstrap interval; the differences between graphs; the same
differences with coverage held constant; Wilcoxon signed-rank tests with
Benjamini-Hochberg correction; and diagnostics on parse failures and evidence
use. Every number in the grounding section of the manuscript comes from the
tables this stage writes.

**Inputs.** `results/tables/G02_run/runs/`,
`results/checkpoints/{G01_benchmark, G02_run}.pkl`, `config/params.yaml` (`grounding`)
**Outputs.** `results/{figures,tables}/G03_analysis/`, `results/checkpoints/G03_analysis.pkl`
**Display item.** Fig G03a (main text), G03b and G03c (supplement)

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | grounding gain per (model, graph); graph contrasts, all pairs and covered-only | parse diagnostics and miss-corrected MRR; rank distribution; `basis`, `evidence_agreement`, `confidence` |

$$\text{gain}_{m,g} = \overline{\text{RR}}_{m,g} - \overline{\text{RR}}_{m,\text{none}}, \qquad \text{contrast}_{m,g,h} = \overline{\text{RR}}_{m,g} - \overline{\text{RR}}_{m,h}$$

where $\overline{\text{RR}}$ is the reciprocal rank averaged over seeds and
shuffles to one value per query, then over queries. A contrast separates when
its 95% bootstrap interval excludes zero and its Wilcoxon *P* survives
BH-FDR within the subset.

## Design notes

### Queries, not rows, are the unit

Nine rows per query and arm are not independent; they share the pair. Rows are
averaged to one reciprocal rank per (model, pair, arm) before anything is
resampled, and the bootstrap resamples pairs. Pairing on `pair_id` (carried on
every row) is exact; the earlier experiment had to replay the sampling to
recover it.

### Two readings of an ordering, and how they are separated

A graph can rank higher because its evidence is better or because it has
evidence more often. Restricting the `main` runs to the `common` pairs holds
availability constant with pools, seeds and models untouched. The shrinkage
in |difference| under that restriction is the share of the apparent ordering
that coverage explains; the residual is what evidence quality is allowed to
claim. The `common` family, run separately, corroborates with an independent
draw.

### Rank test, corrected

Per-query differences are bounded and skewed, so Wilcoxon signed-rank rather
than a paired *t*-test (both are reported). Three graphs give three contrasts
per model and nine per subset; BH-FDR is applied within each subset.

### Failures bias against the best arm

A miss costs more where the arm ranks well, so parse failures make gains
conservative. `mrr_corrected` shows what an arm would score had its failures
ranked like its successes; it is reported, not used.

## Run

    jupyter nbconvert --to notebook --execute grounding_experiment/03_analysis/notebook.ipynb

Seconds to a minute; the bootstrap count is `grounding.n_bootstrap`.

Statistics are `src/biokgsuite/grounding/analysis.py`; figures
`src/biokgsuite/viz/grounding.py`.
