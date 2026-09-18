# 07 · Generalisation

**Computes.** Whether a graph that performs well on average also performs
equitably, on three axes: link prediction restricted to the least-connected
diseases, consistency across therapeutic areas, and prediction of drug-disease
pairs approved after the graph was built. Tier boundaries, area thresholds,
bootstrap count and the random seed come from `config/params.yaml`.

**Inputs.** `config/kgs.yaml`, `config/params.yaml`, `config/kg_snapshots.yaml`,
the KG files named in `kgs.yaml`, `data/gold_standards/bridges/do_diseases.csv`,
`data/gold_standards/bridges/mesh_to_doid.csv`,
`data/gold_standards/bridges/mondo_hierarchy.csv`,
`data/gold_standards/post_cutoff_approvals/post_cutoff_approvals.csv`,
`data/gold_standards/post_cutoff_approvals/post_cutoff_pair_status.csv` (both from
`tools/build_post_cutoff_approvals.py`; see `data/POST_CUTOFF_APPROVALS.md`)
**Outputs.** `results/{figures,tables}/07_generalization/`,
`results/checkpoints/07_generalization.pkl`
**Display item.** Fig 1 panel

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | data-sparse, cross-domain and prospective generalisation | the Q4 → Q1 degradation profile, per-area coefficient of variation, the rank distribution of prospective hits, the network-proximity replication |

$$\text{Generalisation} = \tfrac{1}{3}\left(\text{AUROC}_{Q1} + \overline{\text{AUROC}}_{\text{area}} + \text{AUROC}_{\text{prospective}}\right)$$

A sub-metric whose data is too thin scores N/A and is excluded from the mean,
rather than scoring 0: no signal and no evidence are different claims and a graph
should not be penalised for the second as though it were the first.

## Design notes

### Content dates are evidence, not assumptions

Every prospective claim rests on the graph's content date being right, so each
is now pinned to the publisher's own metadata or to the downloaded artefact, with
the evidence named in `config/kg_snapshots.yaml`. The direction of error matters:
dating a graph too early credits it for an approval it already held, while dating
it too late only discards pairs that were genuinely prospective, so where a range
is defensible the latest date in it is used.

| KG | content date | evidence |
|---|---|---|
| hetionet | 2017-06-08 | v1.0 release record |
| openbilink | 2020-04-27 | Zenodo `publication_date`, 1.0.1 |
| drkg | 2020-05-04 | mtime of `drkg.tsv` inside the published tarball |
| biokg | 2020-06-09 | v1.0.0 release record |
| primekg | 2022-05-02 | Dataverse V2.1 `releaseTime` |
| matrix | 2026-05-28 | Hugging Face commit date |

DRKG is the interesting one: the project publishes no release, tag or checksum,
so the only evidence is the artefact itself, where `drkg.tsv`, `entity2src.tsv`
and `relation_glossary.tsv` all carry 2020-05-04 (the later 2020-06-14 files
under `embed/` are pretrained embeddings, not the graph).

The post-cutoff approval set spans 2023-01-19 to 2026-08-28, so all 128 pairs postdate
every date above except MATRIX's. For the five older graphs the precise day
changes no pair and only has to be established as earlier than 2023-01-19.

### Eligibility is decided before any score is computed

`prospective_max_content_date` in `config/params.yaml` is 2022-12-31. A graph
whose content date falls after it is excluded from this axis and reports N/A
with the reason printed and checkpointed, rather than scoring 0. The rule is
applied to the whole panel in one place, before the evaluation loop, so it
cannot be applied selectively after seeing which graphs it costs.

The threshold follows from the set rather than from any result: the
approvals begin 2023-01-19, so 2022-12-31 is the latest content date at which
every pair in the set is still genuinely prospective. It was set before any
prospective score existed for the four graphs added on 2026-09-18.

NaN and 0 are different claims. A graph that postdates the approvals has not
failed the test; it cannot take it. Averaging a 0 into the dimension score would
make "untestable" and "predicted nothing" indistinguishable, which is the same
defect this benchmark reports in other people's evaluations.

### MATRIX cannot be scored on this axis, and that is the correct result

MATRIX was released on 2026-05-28 and only 5 of the 128 gold pairs postdate it.
Its high raw coverage of the post-cutoff set, 124 of 128 by identifier, reflects
that it already contains these pairs rather than any ability to predict them, and
the date check removes almost all of them. It will fall below `min_new_pairs` and
score N/A, which is what should happen: a graph cannot be evaluated prospectively
against approvals that predate it. Fixing this needs a set of approvals
after 2026-05-28, not a later content date, and the manuscript should say so
rather than leaving a blank cell unexplained.

### Coverage, retrieval and ranking are three questions, not one

A graph can only be scored on the pairs it can express, and the panel differs by
a factor of sixty in how many that is. Worse, the evaluable subsets are not
missing at random: a graph built in 2020 holds the better-studied drugs, which
are also the easier ones, so its conditional score is computed on a friendlier
sample. Comparing two conditional AUROCs across such subsets is comparing
estimates on different populations.

The resolution is to report three numbers and to be exact about which denominator
each one uses.

| | question | denominator |
|---|---|---|
| Coverage | how often can the graph answer? | evaluable / all gold pairs |
| Hits@K (all) | can it deliver a future approval at all? | all gold pairs; inexpressible counts as a miss |
| AUROC (n, CI) | when it answers, is the ranking good? | its own evaluable set |

The asymmetry between the last two is the whole point and it is not arbitrary.
Hits@K is a count whose denominator is the gold pairs, so a pair the graph
cannot rank is a genuine failure to retrieve and needs no imputation. AUROC is
defined over a ranked candidate list, and a pair whose drug has no node has no
list: there is nothing to rank and no neutral value to insert. Imputing zero
would assert that the true positive ranked below every negative, which is a
strong claim rather than an abstention, and it would place the graph beneath
AUROC's 0.5 chance floor, so a graph that cannot represent a drug would score
worse than one guessing at random.

Only the AUROC feeds `temporal_gen` and therefore Fig 1. Coverage is already a
scored dimension (01), and scoring it again here would charge the same
deficiency twice in the composite. Coverage and Hits@K (all) are written to
`prospective_coverage.csv` and carried in the checkpoint so the manuscript can
report them beside the score.

### The pairwise table is the complete-case comparison

The panel-wide intersection is not usable: one low-coverage graph collapses it
to a handful of pairs, and on the current panel it is at most 2. Pairwise
intersections stay large, so each two graphs are compared only on the pairs both
can evaluate, counting how often each ranked the true approval better. Ties are
excluded from the win rate and reported. It assumes nothing about the
distribution of the scores, which is what makes it readable at this sample size,
and a comparison resting on fewer than `min_new_pairs` shared pairs is marked
not comparable rather than quietly reported.

### Disease matching is exact, with no roll-up to a parent term

A prospective pair counts as covered only when the graph holds the disease under
one of the identifiers the set carries. A graph that holds only a
broader parent, or that holds the term with a non-disease type, is not credited.

The case that forced the decision is ischemic stroke (`MONDO:1060198`,
`DOID:0051062`, both valid terms). No graph in the panel holds it: MATRIX has it
as `HP:0002140` typed `biolink:PhenotypicFeature` rather than a disease,
Hetionet has no stroke node at all, and the parent `MONDO:0005098` "stroke
disorder" is present in MATRIX and PrimeKG. Crediting any of those would mean
scoring a graph as able to represent an approval it cannot actually express at
the granularity the approval was granted at, and the roll-up distance would then
be a free parameter chosen after seeing the results. Four of the 101 diseases
are uncovered panel-wide under this rule, and that is reported rather than
repaired.


### Prospective means after the content date and absent from the graph

Two tests, both per graph, both from `tools/build_post_cutoff_approvals.py`
(`data/gold_standards/post_cutoff_approvals/post_cutoff_pair_status.csv`): the approval post-dates
the graph's content date in `config/kg_snapshots.yaml`, and no direct
drug-disease edge of any relation exists between the two nodes. Pairs failing
either are removed from that graph's evaluation and counted, not removed
globally, because a pair can be novel to one graph and asserted in another.
The date test is what keeps a rolling graph such as MATRIX from being scored on
approvals it may already contain; the edge test is what keeps a contemporary
graph from retrieving the answer. Removing the direct edge does not remove the
indirect trace of the approval (co-mention structure, the literature-derived
neighbourhood), so the claim supported is that the answer is not retrievable
and the graph predates the approval, not that the trials that led to it were
invisible.

### A new disease is not a line extension

A PD-1 inhibitor gaining its ninth tumour type is the same mechanism in a
sibling disease. The builder compares each new MONDO term with every disease
the graph already links the drug to, on the MONDO `is_a` hierarchy, and flags
the pair when they are within `prospective_max_hops` (2) steps. The graph's own
indication edges are the reference because they are what the model saw. The
flag is reported by default and becomes a filter with `prospective_strict`.
It can only fire when the graph holds prior indications for the drug, which
`n_prior_indications` in the status table makes visible.

### Coverage is a result, and common support makes graphs comparable

The number of gold pairs a graph can be asked about differs by an order of
magnitude (Hetionet 2, BioKG 66 of 128) because of which drugs and diseases
each graph contains, not because of what it infers. Coverage is written to the
funnel and the checkpoint, and `prospective_common_support` restricts every
graph to the pairs evaluable in all loaded graphs so the AUROC comparison is
on the same questions.

### Intervals are clustered by disease and by drug

Several prospective pairs share a disease, a few drugs contribute many pairs
(trametinib, tislelizumab and upadacitinib 4 each), and each positive carries its own generated negatives, so
rows within a disease or a drug share a candidate pool, a neighbourhood and a
difficulty. A bootstrap that resamples rows treats them as independent and
produces an interval that is too narrow, which is the failure that makes two
graphs look separated when they are not. `bootstrap_*_ci_clustered` resamples
whole clusters instead; on synthetic data with realistic disease-level
difficulty the disease-clustered interval is about 1.4× wider. The row,
disease-clustered and drug-clustered intervals are all recorded in the
checkpoint and the widest is what should be reported.

### Disease resolution goes by identifier, not by name

The post-cutoff approval set is loaded and validated by `post_cutoff_approvals.py`, and diseases
are matched to graph nodes through MONDO, DOID or MeSH, with the name kept only
as a fallback. Name matching was the previous route and it fails silently: a term
the alias map did not anticipate was dropped rather than raised. The notebook
prints how many lookups resolved by identifier against by name, because a set
resolved mostly by name is a weaker set and that should be visible.

### What the loader refuses, and why

A composite `drugbank_id` such as `DB08912,DB08911` matches no graph node, so a
pair carrying one can never be flagged as leaked and never counts as covered: it
passes every filter by matching nothing. That, duplicate (drug, disease) pairs,
a row with no disease identifier in any vocabulary, and unparseable dates are
all fatal.

A disease needs an identifier, but not specifically a DOID. The Disease Ontology
has no term for several real indications in this set (enthesitis-related
juvenile idiopathic arthritis, isolated short stature, complement 3
glomerulopathy) where MONDO does, and MONDO is what PrimeKG keys diseases by, so
requiring DOID would reject rows the pipeline can resolve. 19 of the 128 rows are
MONDO-only and the loader reports that count. A missing MONDO identifier, an approval outside
the declared window, and a DOID whose ontology label shares no word with the
stated indication are reported as warnings, since each is a judgement the reader
should make rather than one the loader should make silently.

### Sparse tiers are defined on total degree

Tiers are quartiles of a disease's total graph degree, so Q1 holds the diseases
with least evidence attached. That is the setting rare conditions sit in, and the
Q4 → Q1 drop is the quantity that says whether a graph's average performance is
carried by well-studied diseases.

### Per-KG random streams

As in 06, every draw is keyed on the graph and the purpose through `_kg_rng`, so
a graph's numbers are the same whether it is run alone or in a panel.

## Run

```bash
# from the repository root, with the six graphs in place (see data/MANIFEST.md)
python tools/run_notebook.py biokgsuite_dimensions/07_generalization/notebook.ipynb
```

Needs all six graphs present: the panel guard refuses a partial run unless
`ALLOW_PARTIAL_PANEL` is set in the config cell. MATRIX dominates the memory
footprint, so run this on the cluster via `hpc/run_main.sbatch` rather than a
laptop. Metric definitions live in `src/biokgsuite/stats/prospective.py` and
`src/biokgsuite/stats/evaluation.py`; the set's loader is
`src/biokgsuite/gold/post_cutoff_approvals.py`, covered by
`tests/test_post_cutoff_approvals.py` and `tests/test_prospective_metrics.py`.

### Original run notes
```bash
# post-cutoff approval set and per-graph checks (once per graph snapshot; loads the graphs)
python tools/build_post_cutoff_approvals.py --extract-views hetionet biokg
python tools/build_post_cutoff_approvals.py --kgs all --figure

jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/07_generalization/notebook.ipynb

python tools/regenerate_figs.py --only 07     # redraw from the checkpoint
python -m pytest tests/test_post_cutoff_approvals.py tests/test_build_post_cutoff_approvals.py tests/test_clustered_bootstrap.py
```

Metric definitions live in `src/biokgsuite/stats/evaluation.py`, the gold-standard
loader in `src/biokgsuite/gold/post_cutoff_approvals.py`, rendering in `src/biokgsuite/viz/dimension07.py`.

## Before the six-graph run

`post_cutoff_pair_status.csv` currently covers only `hetionet` and `biokg`, the
two graphs that existed in `data/` when it was built. Notebook 07 now refuses to
run against a panel it has not all been checked for, because a graph absent from
that table keeps the pairs the content-date and leakage checks would have
removed, making its post-cutoff set larger for a reason that has nothing to do
with the graph.

Rebuild it for the whole panel first:

```bash
python tools/build_post_cutoff_approvals.py --extract-views primekg drkg openbilink matrix
python tools/build_post_cutoff_approvals.py --kgs all
```

The first step loads those four graphs and writes compact views to
`results/checkpoints/07_kg_views/`; MATRIX makes it a cluster job. The second
runs the evaluability checks over every view present and rewrites the status
table. Then run the notebook.
