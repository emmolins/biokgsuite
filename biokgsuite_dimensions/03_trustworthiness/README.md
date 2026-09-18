# 03 · Trustworthiness

**Computes.** How well evidenced and traceable each graph's content is: how
finely an edge can be traced to its origin database, how many independent
upstream resources underpin the graph, and how much of it carries a numeric
confidence.

**Inputs.** `config/kgs.yaml`, `data/SOURCE_COUNTS.md`, the KG files named in
`kgs.yaml`
**Outputs.** `results/{figures,tables}/03_trustworthiness/`,
`results/checkpoints/03_trustworthiness.pkl`
**Display item.** Fig 1 panel

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | edge traceability, source diversity, uncertainty quantification | per-edge provenance examples, Biolink `knowledge_level` breakdown |

$$\text{Trustworthiness} = \tfrac{1}{3}\left(\text{Edge Traceability} + \text{Source Diversity} + \text{Uncertainty Quantification}\right)$$

**Edge traceability** is an ordinal tier read off the graph's own metadata: 0 for
type-level provenance, 0.5 for per-node, 1.0 for per-edge.
**Source diversity** is $\log(n) / \log(n_{\max})$ over the panel, so adding the
tenth source counts for less than adding the second.
**Uncertainty quantification** is the fraction of edges carrying a non-null
numeric confidence field.

## Design notes

### Source counts are documented, not derived

Four of the six graphs carry no per-edge provenance, so no count can be computed
from their edge files. Deriving the number where it is derivable and asserting it
where it is not would make the sub-metric measure metadata format rather than
source breadth. Every figure is instead verified against the resource's own
documentation, recorded with its citation in `data/SOURCE_COUNTS.md`, and
imported from `biokgsuite.stats.provenance.DOCUMENTED` so the number and the
citation cannot drift apart.

`data/SOURCE_COUNTS.md` states the counting rule and lists every source,
including the two counts that required a judgement call (DRKG 7, EC-KG 95).

### Nothing here samples or perturbs the graph

No seed affects the scored metrics. The seed is used only to pick the example
edges shown in the provenance diagnostic, so the scores are exactly reproducible
regardless of `random_seed`.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/03_trustworthiness/notebook.ipynb
```

Metric definitions live in `src/biokgsuite/stats/provenance.py`, rendering in
`src/biokgsuite/viz/dimension03.py`.
