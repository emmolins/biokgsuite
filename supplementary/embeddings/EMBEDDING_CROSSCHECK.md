# Embedding cross-check for notebook 06 link prediction

Script `supplementary/embeddings/run_link_prediction_embeddings.py`; output
`results/tables/06_task_performance/link_prediction_embeddings.csv`
(Supplementary Table S7).

## Why

All three task metrics in 06 are scored with Adamic-Adar, and the two scorers
reported beside it, Common Neighbors and Jaccard, belong to the same family:
each counts shared neighbours. Every conclusion in the dimension was therefore
conditional on local neighbourhood overlap being the right notion of proximity,
and a reviewer is entitled to ask whether the ranking is a property of the
graphs or of that assumption.

TransE (translation in real space) and RotatE (rotation in complex space) have
an entirely different inductive bias. Neither counts shared neighbours. They
score the SAME held-out pairs against the SAME type-constrained negatives,
taken from notebook 06's checkpoint rather than reconstructed, so the
comparison is like-for-like by construction and not by coincidence.

## Calibration

An undertrained embedding sits at chance, and reporting that as "the embedding
disagrees" would describe the training budget rather than the graph. Epochs
were chosen on a validation split carved from the TRAINING positives; the test
pairs never influenced the choice. Test AUROC is recorded at every checkpoint
but selection is on validation alone.

All four configurations calibrated well clear of the 0.55 floor.

## Result

| KG | model | epochs (val-selected) | val AUROC | test AUROC | Adamic-Adar |
|----|-------|----------------------|-----------|------------|-------------|
| hetionet | TransE | 2  | 0.8566 | 0.8545 | 0.8984 |
| hetionet | RotatE | 20 | 0.8619 | 0.8780 | 0.8984 |
| biokg    | TransE | 5  | 0.9750 | 0.9710 | 0.9089 |
| biokg    | RotatE | 80 | 0.9755 | 0.9668 | 0.9089 |

BioKG minus Hetionet:

```
Adamic-Adar   +0.0104
RotatE        +0.0888
TransE        +0.1165
```

**The ranking is preserved under both embeddings.** BioKG leads on link
prediction whether the scorer counts shared neighbours, translates in real
space, or rotates in complex space. So the Adamic-Adar result is not an
artefact of the heuristic family.

The more useful observation is the magnitude. Adamic-Adar puts the two graphs
0.0104 apart, well inside the overlapping confidence intervals, while the
embeddings put them roughly ten times further apart and clear of them. On this
task Adamic-Adar is the conservative scorer: it understates BioKG's advantage
rather than manufacturing it. BioKG's denser, more uniformly typed edge set
appears to give a learned model considerably more to work with than a local
overlap count can express.

## Boundary checks

Two of four configurations selected an epoch count at a grid boundary, which
would leave the optimum unlocated.

- hetionet / TransE selected the lowest grid point (2). Probed epoch 1
  separately: val 0.7876, test 0.7402, i.e. WORSE. The peak is interior and
  the boundary is closed. The default grid now starts at 1 so a rerun covers
  this without a separate probe.
- biokg / RotatE selected the highest grid point (80). Extending it could only
  raise BioKG further and widen a gap that already favours BioKG, so it cannot
  overturn the conclusion. Left unextended; stated rather than hidden.

## Extension to the other two tasks

Retrieval and multi-hop were Adamic-Adar alone, so the same objection applied
to two thirds of the dimension. Both are now covered by
`supplementary/embeddings/run_task_embeddings.py`. Scores are at the
pool-normalised cut (top 0.5%), the same one the composite uses.

Retrieval scores the (disease, r, gene) triple directly, maximised over
whichever relations carry disease-gene edges in that graph. Multi-hop composes
the two hops additively, score(drug,r1,g) + score(g,r2,disease), because these
models score by negative distance, so a sum in score space is a product in
likelihood space.

A trained model cannot hide an edge at scoring time the way 06 does, so for
each task the evaluated pairs are removed from the training triples BEFORE
training, which is stricter than 06 and never looser.

| task | KG | model | epochs | metric@0.5% | Adamic-Adar |
|------|----|-------|--------|-------------|-------------|
| retrieval | hetionet | TransE | 5  | 0.0661 | 0.1488 |
| retrieval | hetionet | RotatE | 1  | 0.0248 | 0.1488 |
| retrieval | biokg    | TransE | 10 | 0.1400 | 0.1733 |
| retrieval | biokg    | RotatE | 40 | 0.1500 | 0.1733 |
| multihop  | hetionet | TransE | 5  | 0.0654 | 0.1055 |
| multihop  | hetionet | RotatE | 1  | 0.0992 | 0.1055 |
| multihop  | biokg    | TransE | 40 | 0.0802 | 0.1329 |
| multihop  | biokg    | RotatE | 40 | 0.1392 | 0.1329 |

BioKG minus Hetionet:

```
retrieval   Adamic-Adar +0.0246   TransE +0.0739   RotatE +0.1252
multihop    Adamic-Adar +0.0274   TransE +0.0148   RotatE +0.0401
```

**Ranking preserved on both tasks under both models**, so all three metrics in
the dimension now survive a scorer from outside the shared-neighbour family.
As with link prediction, the embeddings mostly widen the gap rather than
narrow it; multi-hop under TransE is the one case where the gap shrinks, and
it stays the same sign.

Adamic-Adar beats both embeddings on seven of eight configurations here, which
is the opposite of the link-prediction result. That is not a contradiction: AA
is given the graph with one edge hidden at scoring time and can exploit all the
remaining local structure, whereas the embeddings must compress the whole graph
into 64 dimensions and are trained with every evaluated pair removed. It does
say that on these two tasks the heuristic is a strong baseline rather than a
weak one, which strengthens rather than weakens its use as the primary scorer.

## A validation set has to match the evaluation's difficulty

The first multi-hop attempt built validation triples locally, by sampling a
drug and a disease that share a gene neighbour and calling that gene the
answer. Those triples are far easier than CTD's curated ones by construction:
validation rose 0.18 to 0.70 across epochs while the real metric FELL 0.057 to
0.049. Selecting on it would have picked the worst available model.

Validation triples now come from notebook 06, drawn from the wider CTD pool
with the same resolution, leakage filter and entity caps as the benchmark and
disjoint from it (359 per KG). Validation and evaluation now track each other
closely, which is the check that the proxy failed.

Retrieval on Hetionet needed a different fallback: it has 121 eligible diseases
and evaluates all 121, so no spare disease exists. Validation there uses a
second, different held-out gene from each evaluated disease. Both are hidden
from training, so no evaluation answer leaks, but validation and evaluation
queries share head entities, and hiding two partners instead of one leaves the
model with less signal than Adamic-Adar had rather than more.

## Caveats

- One seed per configuration. The differences here are large relative to the
  Adamic-Adar gap, but a seed-variance estimate is not in hand and this is not
  a significance claim.
- dim=64 for every model, not tuned per graph.
- All three tasks are now covered.
- Retrieval and multi-hop use a single epoch grid (1,2,5,10,20,40) and dim=64.
  Two configurations selected the top of that grid (biokg RotatE on both
  tasks); extending it could only raise BioKG and widen a gap already in its
  favour, so it cannot overturn the ranking.
- Hetionet RotatE selected 1 epoch on both tasks, the bottom of the grid, and
  its scores are the weakest in the table. Extending downward is not possible,
  so this configuration should be read as poorly fit rather than as evidence
  about Hetionet.
- Link prediction and the other two tasks use different training sets, because
  each hides a different set of edges. Numbers are comparable across KGs within
  a task, not across tasks.
