# 01 · Coverage

**Computes.** How completely each KG represents the entities and relations that
drug repurposing needs, against task-derived reference sets that none of the six
graphs ingests.

**Inputs.** `config/kgs.yaml`, `config/params.yaml`, `data/reference_sets/`,
the KG files named in `kgs.yaml`
**Outputs.** `results/{figures,tables}/01_coverage/`,
`results/checkpoints/01_coverage.pkl`
**Display item.** Fig 1 panel; Table S6

## What is scored

| | scored | computed and stored, not scored |
|---|---|---|
| Entities | Drug, Disease, Gene/Protein | Pathway |
| Relations | drug→gene, disease→gene, disease→disease, drug→disease | drug→pathway, disease→pathway |

$$\text{Coverage} = \tfrac{1}{2}\left(\text{Entity Coverage} + \text{Relation Coverage}\right)$$

Each sub-metric averages over scored types that are evaluable in **every** KG,
so the comparison is like for like.

Disease coverage is additionally reported by stratum (all / rare / with an
approved indication / without). That split is the headline result, not a
robustness check: it is what shows whether the graphs are populated where
repurposing is most needed.

## Design notes

### Pathways are excluded for vocabulary mismatch, not irrelevance
 The six graphs use Reactome, KEGG, SMPDB and PathwayCommons
accessions, so matching on identifier asks which source a graph ingested rather
than what biology it holds. The notebook's pathway diagnostic prints the
evidence: Hetionet has 1,822 pathway nodes, 0 carrying a Reactome accession, yet
873 of its 1,808 distinct pathway names (48.3%) are exact matches to Reactome
pathway names. That 48.3% is a floor, since the reference holds only the
Reactome pathways with a human drug target; a one-off check against
PathwayCommons v7 put the true agreement near three quarters. Full argument, the
PC7 check, and the route back via gene-set matching are in
`data/REFERENCE_SETS.md` under "Why pathways are not scored".

Nothing is discarded. The pathway reference sets are still built and
version-stamped, and the checkpoint carries `pathway_diagnostic`,
`scored_entity_types`, `scored_relations` and `reported_not_scored`. Re-including
pathways means editing `SCORED_*` in the notebook, not rebuilding data.

### Resolvability vs coverage

Coverage alone conflates "the KG does not contain this" with "we could not map
its identifier". The notebook prints a resolvability table first: the fraction
of each KG's own IDs that reach the reference namespace. A low value there caps
the coverage that graph can attain for pipeline reasons rather than content
reasons, and the relation gap decomposition splits every missing gold pair into
entity-absent, unresolved, or relation-absent.

## Run
```bash
jupyter nbconvert --execute --inplace notebook.ipynb
# or, headless:
python tools/run_notebook.py biokgsuite_dimensions/01_coverage/notebook.ipynb
```

Reference sets must exist first:
```bash
bash   tools/fetch_opentargets.sh      # Open Targets 26.06
python tools/build_crosswalks.py       # MONDO v2026-08-04, DOID v2026-08-31
python tools/build_reference_sets.py   # writes data/reference_sets/
```

Metric definitions live in `src/biokgsuite/stats/coverage.py`, rendering in
`src/biokgsuite/viz/primitives.py`.
