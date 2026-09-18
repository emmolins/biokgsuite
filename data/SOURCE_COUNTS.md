# Source counts for the trustworthiness dimension

Verified 2026-09-15. These feed `source_diversity` in notebook 03 and are
hardcoded in `src/biokgsuite/stats/provenance.py` (`DOCUMENTED`).

## Counting rule

A **source** is a distinct upstream body of evidence that contributes edges, as
enumerated by the resource itself.

- **Literature-derived corpora count.** They are independent bodies of evidence
  with their own failure modes. Excluding them would also be inconsistent with
  EC-KG's 95, which reaches SemMedDB (literature-mined) through RTX-KG2.
- **Ontologies and identifier-mapping resources do not count** when they supply
  vocabulary rather than edges. This is why OpenBioLink is 9 and not 16.
- **Aggregators are counted as the de-duplicated union of their primary
  sources** where the resource publishes that union, and as themselves where it
  does not.

The rule cannot be applied at one level for every graph, because the six sit at
different depths of the same aggregation hierarchy: DRKG names Hetionet as one
of its sources, and EC-KG consumes PrimeKG. Where a release carries per-edge
provenance, `provenance.observed_sources` derives the count from the data and
notebook 03 reports it beside the documented figure. Where the two disagree,
the disagreement is the finding.

## Counts

| KG | n | Basis | Citation |
|---|---|---|---|
| Hetionet | 29 | "We obtained and integrated data from 29 publicly available resources to create Hetionet v1.0" | Himmelstein et al. 2017, *eLife* 6:e26726. https://elifesciences.org/articles/26726 |
| PrimeKG | 20 | "PrimeKG integrates 20 high-quality resources to describe 17,080 diseases" | Chandak, Huang & Zitnik 2023, *Sci Data* 10:67. https://www.nature.com/articles/s41597-023-01960-3 |
| DRKG | **7** | 6 databases + Bibliography (see note) | DRKG README statistics tables. https://github.com/gnn4dr/DRKG |
| BioKG | 13 | 13 of the 18 databases in Table 1 are ticked in the "In BioKG?" column | Walsh, Mohamed & Nováček 2020, *CIKM '20*. https://dl.acm.org/doi/10.1145/3340531.3412776 |
| OpenBioLink | 9 | 9 distinct edge-contributing sources in the README source table | https://github.com/OpenBioLink/OpenBioLink ; Breit et al. 2020, *Bioinformatics* 36:4097 |
| MATRIX / EC-KG | **95** | "95 aggregated knowledge sources" | Kaniewski et al. 2026, *bioRxiv* doi:10.64898/2026.08.26.747253 |

## Notes on the two contested counts

### DRKG is 7, not 6

The README prose reads: "DRKG includes information from six existing databases
including DrugBank, Hetionet, GNBR, String, IntAct and DGIdb, **and** data
collected from recent publications particularly related to Covid19."

That sentence puts the literature outside the six because a corpus is not a
database. The README's own statistics tables do not:

| | Drugbank | GNBR | Hetionet | STRING | IntAct | DGIdb | **Bibliography** | Total |
|---|---|---|---|---|---|---|---|---|
| Entities | 19,911 | 44,033 | 45,279 | 18,316 | 16,474 | 8,899 | **9,464** | 97,238 |
| Triplets | 1,424,790 | 335,369 | 2,250,197 | 1,496,708 | 256,151 | 26,290 | **84,756** | 5,874,261 |

Bibliography has its own column and contributes 84,756 edges, including 58,629
gene-gene and 25,666 compound-gene. Under the counting rule above it is a
seventh source.

Kaniewski et al. 2026 report 6 for DRKG, following the prose rather than the
tables. Cite ours as 7 and note the discrepancy rather than hiding it.

### MATRIX is EC-KG, and the paper supersedes the docs page

Every Cure's tech docs list five first-level knowledge sources (ROBOKOP,
RTX-KG2, SPOKE, EmBiology, PrimeKG). The paper describes the released graph as
integrating **three** (RTX-KG2, ROBOKOP, PrimeKG) and names SPOKE as future
work: "The integration framework is extensible and can integrate additional
knowledge graphs, including large-scale resources such as SPOKE."

95 is the de-duplicated union of those three graphs' primary sources, not a sum:
"only 15 out of 95 aggregated sources (16%) overlap between RTX-KG2 and ROBOKOP
KG, with PrimeKG being a distinct entity with no formal overlap."

The graph should be called EC-KG in the manuscript; MATRIX is the platform.

## Cross-check

Kaniewski et al. 2026 independently report counts for three of our panel, all
matching the primary publications, which is a useful check on all three:

| KG | their n | nodes | edges |
|---|---|---|---|
| EC-KG | 95 | 7,354,612 | 81,149,812 |
| Hetionet | 29 | 47,031 | 2,250,197 |
| PrimeKG | 20 | 129,375 | 4,050,249 |
| DRKG | 6 | 97,238 | 5,874,261 |

Their Hetionet node and edge counts also match what `load_kg` reads from our
copy of the release.
