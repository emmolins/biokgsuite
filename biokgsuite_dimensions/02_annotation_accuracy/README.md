# 02 · Annotation accuracy

**Computes.** Whether each graph's own labelling survives a check from outside
the graph. The graph asserts a node is a Drug, or that an edge joins two Drugs;
`IdResolver` is asked whether the identifier bears that out.

**Inputs.** `config/kgs.yaml`, `config/id_schemes.yaml`, `data/crosswalks/`,
the KG files named in `kgs.yaml`
**Outputs.** `results/{figures,tables}/02_annotation_accuracy/`,
`results/checkpoints/02_annotation_accuracy.pkl`
**Display item.** Figs 2 and 3 (via `00_summary`); Supplementary Table S3a

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | entity type accuracy, relation endpoint accuracy | node index uniqueness, relation exclusivity |

$$\text{Annotation accuracy} = \tfrac{1}{2}(\text{entity type accuracy} + \text{relation endpoint accuracy})$$

**Entity type accuracy.** Of the nodes whose identifier resolves at all, the
fraction resolving into the namespace matching their declared type. A node
labelled `Drug` whose identifier resolves only into MONDO is namespace leakage.

**Relation endpoint accuracy.** Of the edges whose endpoints both resolve, the
fraction where both resolve into the types the edge declares, so a `drug_drug`
edge really does join two drugs.

Relation exclusivity stays a diagnostic because a graph declaring no mutually
exclusive relations would take a vacuous 100%, rewarding schema minimalism over
curation.

## Design notes

### Why the arbiter is external

The previous implementation checked each graph against a schema hand-written by
reading that graph, and DRKG's schema was parsed out of its own relation names
while MATRIX's accepted every type pair. Three of six graphs could not fail.
Worse, the denominator counted edges whose relation was absent from the schema,
so an incomplete schema *raised* the score: the less schema work a graph
received, the better it looked.

Nothing is now asserted about what a relation ought to join. The edge's own
`x_type` and `y_type` are the claim; the identifiers are asked whether the claim
holds; every graph is held to the same rule.

### Reading the numbers

Report `Accuracy (%)` next to `nodes_checkable_pct` and `edges_checkable_pct`.
A perfect score over a quarter of the edges is a weaker claim than the same
number over most of them. Edges are excluded when an endpoint is typed outside
the four canonical types (side effects, anatomy, GO terms) or fails to resolve.
Hetionet, for example, scores 1.00 over 86% of nodes but only 24.9% of edges.

Unresolvable identifiers are excluded from both denominators rather than
counted as wrong. That is a gap in our crosswalks, which nb01 reports as
resolvability, not a mislabelling by the curator.

### Bare numbers

Check 1 reads the shape of an identifier, and a bare number has almost none:
`8150` fits MONDO, Entrez, GO, UBERON and HPO alike. PrimeKG writes five of
its node types that way (its `x_source` column says which registry each is:
MONDO for diseases, NCBI for genes, GO, UBERON, HPO), Hetionet and DRKG write
genes as `Gene::1956`, and the MATRIX loader emits bare MONDO and Entrez
numbers after bridging. `config/id_schemes.yaml` therefore declares the bare
shape under every type that a graph writes it for, marked ambiguous, so the
declared type decides how the number is read, and the reading is then checked
against the registry rather than assumed. `tools/audit_identifiers.py`
reports, per graph and declared type, every identifier shape and the share
each registry confirms. On PrimeKG the bare numbers confirm at 94.6%
(MONDO), 98.8% (Entrez), 84.5%, 90.0% and 98.9% (GO biological process,
molecular function, cellular component), 99.8% (UBERON) and 100% (HPO). The
shortfall is retired terms: PrimeKG was built on 2022 releases, and a term
retired since resolves to its `replaced_by` successor when one is named and
is otherwise unconfirmed. Before the bare shapes were declared, every one of
those nodes was judged a mislabelled gene and 47% of PrimeKG's edges were
counted as type violations, which put its annotation accuracy at 0.76.

### Live API spot-check

Removed from the notebook. It hit eight external services, cached transient
failures permanently as negative results, and slept on cache hits, while
checking the same identifiers the resolver now checks offline against dated,
checksummed snapshots. If live confirmation is wanted, it belongs in a
standalone tool whose output is a dated table, not in the scored path.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/02_annotation_accuracy/notebook.ipynb
```

Metric definitions live in `src/biokgsuite/stats/annotation.py`, rendering in
`src/biokgsuite/viz/primitives.py`.
