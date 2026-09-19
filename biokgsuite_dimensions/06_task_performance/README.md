# 06 · Task performance

**Computes.** Whether a graph's structure supports biological inference, on three
tasks run with neighbourhood heuristics only: link prediction on drug-disease
indication pairs, neighbourhood retrieval of held-out disease-gene associations,
and multi-hop recovery of mechanistic target genes. No embeddings and no
training, so what is measured is structural signal rather than a model fitted to
the graph. Sample sizes, bootstrap count and the random seed all come from
`config/params.yaml`.

**Inputs.** `config/kgs.yaml`, `config/params.yaml`, the KG files named in
`kgs.yaml`, `data/gold_standards/ctd/CTD_chem_gene_ixns.tsv.gz`,
`data/gold_standards/ctd/CTD_curated_genes_diseases.tsv.gz`
**Outputs.** `results/{figures,tables}/06_task_performance/`,
`results/checkpoints/06_task_performance.pkl`
**Display item.** Fig 1 panel

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | link prediction, neighbourhood retrieval, multi-hop reasoning | Common Neighbors and Jaccard, conditional Hits@100, the embedding cross-check |

$$\text{Task performance} = \tfrac{1}{3}\left(\frac{\text{AUROC} - 0.5}{0.5} + \text{Recall@100} + \text{Hits@100}\right)$$

Link prediction is AUROC under Adamic-Adar with type-constrained negatives.
Neighbourhood retrieval is Recall@100 from leave-one-out over diseases with at
least five known targets. Multi-hop is Hits@100 over known (drug, gene, disease)
triples scored by AA(drug, *g*) × AA(*g*, disease).

## Design notes

### The three sub-metrics do not share a chance floor

AUROC is 0.5 under the null; Recall@100 and Hits@100 are 0. Averaging the raw
values lets a graph with no link-prediction signal contribute a third of a point
to the dimension, which is not a claim any of the three measurements support.
Link prediction is therefore mapped to (AUROC − 0.5) / 0.5, which is 2·AUROC − 1,
the Gini coefficient, before the mean,
under the `RESCALE_AUROC` flag, so all three run from 0 at chance to 1 at
perfect. The unscaled AUROC is still what the figure and the table report,
because that is the number comparable to the literature.

### Adamic-Adar is pre-specified, and the other two are the check on it

Adamic-Adar down-weights shared neighbours that are high-degree hubs
($\sum_{w \in N(u) \cap N(v)} 1/\log|N(w)|$), so a low-degree shared intermediate
counts as more specific evidence than a hub. That is a biological argument, which
is why it is chosen in advance rather than after seeing the results. Common
Neighbors and Jaccard run alongside: if the three move together the conclusion is
a property of the graph, and if Adamic-Adar stands apart it is a property of its
weighting, which the manuscript then has to say.

### Conditional Hits@100 separates two different failures

A multi-hop triple can miss because the intermediate gene ranks badly, or because
it carries no score at all. Only the first is a ranking failure; the second is
absence of evidence. The unconditional figure merges them, so the conditional
form, computed over triples that received any score, is reported beside it.

### Hits@K saturates on a small candidate pool

With `NEG_RATIO` negatives per positive the scored set is
(1 + NEG_RATIO) × n_sampled, so Hits@100 approaches 1 by construction once that
set is near 100. The notebook prints a warning when the pool is that small, and
the rank denominator recorded in the checkpoint is the scored candidate set, not
the number of diseases in the graph.

### Per-KG random streams

Every draw is keyed on the graph and on what the draw is for, through `_kg_rng`.
A single generator consumed across the KG loop makes each graph's splits,
negatives and intervals depend on how many graphs preceded it, so adding a graph
to `config/kgs.yaml` would silently move the numbers of every graph after it.

### The embedding cross-check

The same three tasks are rerun under a learned representation. Agreement with the
heuristic ranking means the conclusion is structural rather than an artefact of
one scorer; disagreement is itself a result, and either way it is reported rather
than folded into the score.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/06_task_performance/notebook.ipynb

python tools/regenerate_figs.py --only 06    # redraw from the checkpoint
```

Metric definitions live in `src/biokgsuite/stats/evaluation.py` and
`src/biokgsuite/graph/scoring.py`, rendering in
`src/biokgsuite/viz/dimension06.py`.
