# 05 · Stability

**Computes.** Robustness of Adamic-Adar rankings of known drug-disease pairs to
edge dropout, random and periphery-first, against a degree-matched null. Sample
size, dropout rates, trial count, top-K window and the primary scoring rate all
come from `config/params.yaml: analysis_params.robustness`.

**Inputs.** `config/kgs.yaml`, `config/params.yaml`, the KG files named in
`kgs.yaml`
**Outputs.** `results/{figures,tables}/05_stability/`,
`results/checkpoints/05_stability.pkl`
**Display item.** Fig 1 panel

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | random-dropout stability, periphery-dropout stability | the degree-matched null, Common Neighbors and Jaccard, top-K overlap |

$$\text{Stability} = \tfrac{1}{2}\left(\text{Random-dropout stability} + \text{Periphery-dropout stability}\right)$$

Each sub-metric is $(r + 1) / 2$, a linear rescaling of Spearman *r* on therapeutic pairs
under Adamic-Adar at `primary_dropout`. The full degradation profile across every
configured rate is written to `results/tables/05_stability/stability.csv` and
drawn in the figure.

**Thresholds** (interpretive, not used in scoring): *r* > 0.95 at the primary
rate = high stability; *r* < 0.85 = fragile topology.

## Design notes

### The null is what carries the biological claim

Adamic-Adar is a function of shared-neighbour degrees, so a high correlation on
therapeutic pairs alone does not establish that a graph encodes robust
*therapeutic* structure: it may only show that hub-driven rankings survive edge
loss. The null repeats the whole measurement on random drug-disease pairs whose
endpoint degrees match the therapeutic pairs, and the gold-minus-null gap is what
any claim about biological signal has to be read against.

The null is reported rather than scored, in the same way pathway evidence is
reported in notebook 01 and NMI/ARI in notebook 04; the sub-metrics remain
the absolute rescaled correlation on the therapeutic pairs. Null correlations rest on fewer pairs
than gold, because degree-matched random pairs more often share no neighbours at
all; both counts are recorded in `stability.csv`.

### Both strategies are replicated

Random dropout is stochastic by construction. Periphery dropout is deterministic
given a tie-break, but two edges with the same minimum endpoint degree are
equally peripheral by this metric's own definition, so which goes first is
arbitrary rather than a fact about the graph. Each trial shuffles inside the
degree band under a seeded permutation, so the reported SD says how much the
periphery score depends on the draw. Every sub-metric therefore carries an
uncertainty; previously half the score had none.

### Periphery order is fully determined by the graph

Edges are sorted by (minimum endpoint degree, sum of endpoint degrees, sorted
endpoint indices). Sorting on minimum degree alone leaves tie groups large enough
that most of what is removed would be decided by the order the loader built the
graph in: at the primary rate on BioKG the boundary degree is 4 and its tie group
holds tens of thousands of edges.

### Dropout is spent on assertions, not on node pairs

`build_graph` produces a simple graph, because Adamic-Adar, clustering,
modularity and LCC are all defined that way, but each edge records the full set
of predicates asserting it and how many triples it carries. Dropout draws
uniformly from the triples and an edge disappears only when every triple
supporting it is gone. Removing whole pairs instead made a triply-asserted pair
exactly as easy to sever as a singly-asserted one, severed with probability *r*
rather than roughly *r*³, which does not match the incompleteness the dimension
claims to model. The periphery budget is counted in triples too, so both
strategies remove the same amount of evidence at the same nominal rate.

One consequence is worth reading off the output rather than assuming: an equal
triple budget severs an *unequal* number of node pairs under the two strategies,
because redundantly asserted pairs are not distributed evenly across the degree
range. The sweep prints both counts at every rate.

### Three scorers, one scored

The sweep runs under all three neighbourhood heuristics in
`biokgsuite.graph.scoring`: Adamic-Adar (scored, as Table 1 states), Common
Neighbors and Jaccard (reported). If the three move together, the conclusion is a
property of the graph; if Adamic-Adar stands apart, it is a property of its
inverse-log-degree weighting, and the manuscript has to say so.

Top-K overlap (K = `robustness.top_k`) is reported for characterisation and does
not enter the score. It is undefined, and reported as such, when K covers the
whole pair sample.

### Reading the numbers

Reported sample sizes are the pairs actually correlated, not the pairs sampled.
Pairs with no Adamic-Adar signal before or after perturbation carry no rank
information and are excluded; `sample_size` pairs are drawn but fewer are
correlated, and the stability table records both.

The rescaling compresses, so the scale is narrow at the top and this dimension separates
graphs far less than the others do. The degradation profile across rates is more
informative than the single scored number.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/05_stability/notebook.ipynb

python tools/regenerate_figs.py --only 05    # redraw from the checkpoint
python -m pytest tests/test_stability.py     # metric unit tests
```

Metric definitions live in `src/biokgsuite/stats/stability.py`, rendering in
`src/biokgsuite/viz/dimension05.py`.
