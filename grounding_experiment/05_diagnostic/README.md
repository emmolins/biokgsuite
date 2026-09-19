# G05 · Diagnostic

**Computes.** Per-pair difficulty on the KG arms, whether the models agree on
which pairs are hard, and what predicts difficulty: therapeutic area, how many
graphs cover the pair, the model's reported confidence and its reported
agreement with the evidence.

**Inputs.** `results/tables/G02_run/runs/` (the headline models, `main`
family), `results/checkpoints/G01_benchmark.pkl`, `config/params.yaml` (`grounding`)
**Outputs.** `results/{figures,tables}/G05_diagnostic/`, `results/checkpoints/G05_diagnostic.pkl`
**Display item.** Supplementary figure (hard pairs)

## What is measured

| | measured | reported, not scored |
|---|---|---|
| | per-pair MRR (all KG arms, grounded arms only, `none`); cross-model correlation | tertile summaries, covariate correlations, area means |

Nothing here is a score; the stage explains the G03 result.

## Design notes

### Difficulty belongs to the pair if the models agree

If per-pair MRR correlates strongly across models, the hard pairs are hard
for a reason in the pair (an indirect mechanism, a disease the graphs know
little about), not in the model. That is what the cross-model matrix tests.

### Coverage attributed fairly

`mrr_grounded` averages only the arms whose graph had evidence for the pair,
so a pair looks hard only when the evidence was there and did not help, not
because two graphs had nothing to say. The number of graphs with evidence is
the covariate that predicts difficulty in the earlier experiment; the figure
keeps both readings side by side, and `n_kg_covered` stays in the table so
the two counts can be compared.

### Covariates come from the benchmark, not from names

Modality (Open Targets drug type), area (MONDO hierarchy) and approval year
are columns of the benchmark frame; nothing is inferred from a drug's name
inside this stage.

## Run

    jupyter nbconvert --to notebook --execute grounding_experiment/05_diagnostic/notebook.ipynb

`difficulty_table` and `cross_model_corr` are in `src/biokgsuite/grounding/analysis.py`;
figures in `src/biokgsuite/viz/grounding.py`.
