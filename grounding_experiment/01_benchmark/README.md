# G01 · Benchmark

**Computes.** The benchmark every later stage reads: the 128 post-cutoff
approval pairs with their per-graph coverage, leak and evaluability flags, the
`common` subset evaluable in every KG arm, the candidate pools for the first
seed, and a selection audit asking whether `common` keeps a different kind of
pair from the rest.

**Inputs.** `config/params.yaml` (`grounding`),
`data/gold_standards/post_cutoff_approvals/{post_cutoff_approvals.csv, post_cutoff_pair_status.csv}`,
`data/gold_standards/bridges/{approved_drug_pool.csv, approved_indications.csv, mondo_hierarchy.csv, do_diseases.csv, mesh_to_doid.csv}`
**Outputs.** `results/{figures,tables}/G01_benchmark/`, `results/checkpoints/G01_benchmark.pkl`
**Display item.** Benchmark overview (`G01d_benchmark_overview`: approvals
against each graph's content date, per-graph evaluability, how many graphs can
evaluate each pair), intended for the main text; `G01a` to `G01c` (selection
audit) for the supplement

## What is measured

| | measured | reported, not scored |
|---|---|---|
| | per-graph covered / grounded / leaked / evaluable counts; size of the `common` subset | benchmark composition, distractor-pool statistics, the selection audit |

Nothing here is a score. `covered_<kg>` is the builder's flag: the graph
holds the drug at the degree floor and the disease. `grounded_<kg>` is the
builder's second flag, computed by the dossier code itself: the drug has at
least one gene or pathway neighbour a dossier would render, so it does not
reach the model blank. `evaluable_<kg>` needs both, plus an approval that
postdates the graph's content date and no direct edge; with `grounding.strict`
it also excludes pairs within `prospective_max_hops` of a prior indication.
Notebook 07 uses the same status table without the grounded flag, because a
link-prediction scorer needs the node, not evidence text; the grounding
experiment folds it in because a blank dossier is what the model is scored on.

### Covered is not grounded

A drug can clear the degree floor on drug-drug interaction or side-effect
edges alone. Such a drug is in the graph, and notebook 07 can score it, but a
dossier (gene and pathway neighbours only, ADME relations dropped where the
graph labels them) has nothing to render, and the positive reaches the model
as a bare label beside distractors that carry evidence. The two flags are kept
apart, the gap between them is tabulated per graph (`graph_status.csv`), and
G03 measures what that gap costs (the absence penalty) instead of folding it
into the pooled lift unseen.

## Design notes

### One benchmark for both parts of the paper

The positives are exactly the pairs notebook 07 scores, and the per-graph
flags are read from the status table `tools/build_post_cutoff_approvals.py`
writes rather than recomputed, so the coverage-controlled comparison in G03
uses the same checks the generalisation dimension reports. The source table's
SHA-256 travels in the checkpoint.

### Two query families from the same pairs

`main` asks every pair in every arm and controls for coverage afterwards by
restricting to the `common` pairs. `common` asks only those pairs, as its own
run with its own distractor draw. The first is the within-experiment control;
the second corroborates it with an independent sample. Both are cheap enough
to run because `common` is a subset.

### Distractors that are wrong under either record

A distractor is an approved drug (Open Targets, DrugBank id via cross-reference
or the name bridge, `tools/build_distractor_pool.py`) that is not indicated for
the target disease or any MONDO descendant of it, and is not another
post-cutoff positive for that disease. The descendant closure is new: a drug
approved for a subtype is a defensible answer for the parent term and must not
be a distractor. Distractors are drawn without regard to graph membership; the
run records how many candidates carried evidence so G03 can check that a
positive was not simply the only candidate with a dossier.

### Order is seeded and stratified

Queries are permuted within therapeutic area and taken round-robin across
areas, so a run cut short by a rate limit is not an oncology run. The area
comes from the MONDO hierarchy (`OntologyClassifier.classify_mondo`), with the
MeSH C-category terms on the name as the fallback; 14 pairs remain
unclassified and form their own bucket.

### The selection audit

`kept` marks the `common` pairs. Fisher odds ratios per categorical level and
Mann-Whitney tests per continuous covariate say whether evaluability in every
graph tracks modality, area, year or how many regulators approved the pair.
A trait that predicts being kept narrows what the coverage-controlled result
generalises to; the audit is reported so a reader can weigh that.

## Run

    jupyter nbconvert --to notebook --execute grounding_experiment/01_benchmark/notebook.ipynb

Runs in seconds; no graph is loaded.

Benchmark construction lives in `src/biokgsuite/grounding/benchmark.py`, the
audit statistics in `src/biokgsuite/grounding/bias.py`, the figures in
`src/biokgsuite/viz/grounding.py`.
