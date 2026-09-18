# 04 · Topology

**Computes.** Whether each graph's structure reflects meaningful biological
organisation: how much of it hangs together, how locally clustered it is against
a degree-matched random baseline, whether validated drug-disease pairs survive
edge removal better than matched random pairs, and how well its communities line
up with entity types.

**Inputs.** `config/kgs.yaml`, `config/params.yaml`, `data/reference_sets/`,
`data/crosswalks/`, the KG files named in `kgs.yaml`
**Outputs.** `results/{figures,tables}/04_topology/`,
`results/checkpoints/04_topology.pkl`
**Display item.** Fig 1 panel

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | graph cohesion, clustering coefficient, graph reachability, community purity | NMI and ARI beside AMI, composition of the largest communities |

$$\text{Topology} = \tfrac{1}{4}\left(\text{Cohesion} + \text{Clustering} + \text{Reachability} + \text{Purity}\right)$$

## Design notes

### Community purity is scored with AMI, not NMI

AMI corrects for chance. NMI rises with the number of communities even when the
partition carries no information, which would reward a graph for being
fragmented. NMI and ARI are computed and stored beside AMI as diagnostics.

### AMI is graded against canonical entity types, not each graph's raw labels

BioKG splits diseases into `Disease` (MeSH) and `GeneticDisorder` (OMIM);
`config/kgs.yaml` maps both to `Disease` and notebook 01 counts them as one type.
Grading on the raw split would penalise a graph for a labelling artefact this
benchmark has already adjudicated, and would let the graph's own schema decide
what counts as a pure community, which is the arbiter problem this framework
exists to avoid. A raw label the config does not map keeps its own class and is
named in the output rather than folded in silently.

### Clustering is sampled for every graph, at one fixed size

The cost of the exact computation is set by the degree distribution, not the node
count: it sums d(v)² over every node, so it is dominated by hubs and does not
scale to the larger graphs in the panel. A fixed k-node sample replaces that sum
with k draws from the same distribution, so its cost is bounded by the sample
size regardless of how large the graph is. Applying one sample size to every
graph also keeps the comparison like for like, rather than computing some exactly
and estimating others. The standard error the sample induces on the score is
reported per graph.

### Reachability is differential, and shares notebook 01's gold standard

The gold standard is `data/reference_sets/relation_drug_disease.tsv` (Open
Targets approved indications), resolved through the shared `IdResolver`. The
score is the gap against type-matched random pairs, so a larger or denser graph
that preserves *all* pairs equally gains nothing unless the gold pairs are
specifically more resilient.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/04_topology/notebook.ipynb
```

Metric definitions live in `src/biokgsuite/stats/topology.py`, rendering in
`src/biokgsuite/viz/dimension04.py`.
