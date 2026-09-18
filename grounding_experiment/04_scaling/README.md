# G04 · Scaling

**Computes.** The grounding gain at each rung of a same-family model ladder
(Llama 3.x, 1B to 405B), on the `main` family, split into the lift pooled over
every pair and the lift on the pairs the graph actually covers.

**Inputs.** `results/tables/G02_run/runs/` (the models in
`grounding.capacity_models`), `results/checkpoints/G01_benchmark.pkl`,
`config/params.yaml` (`grounding`)
**Outputs.** `results/{figures,tables}/G04_scaling/`, `results/checkpoints/G04_scaling.pkl`
**Display item.** Supplementary figure (capacity series)

## What is scored

| | scored | reported, not scored |
|---|---|---|
| | MRR per arm per rung; pooled and covered lift per rung | the per-graph breakdown by rung |

## Design notes

### Why a within-family ladder

The cross-family comparison in G03 confounds capability with architecture and
training data. A ladder within one family holds those roughly fixed. No Llama
generation offers five sizes, so the rungs span generations (1B and 3B from
3.2, 8B and 405B from 3.1, 70B from 3.3); the generation is a column in the
table so the confound stays visible.

### Pooled versus covered

A larger model has memorised more approvals, which raises the `none` arm and
shrinks the pooled lift mechanically. A graph that does not contain the pair
cannot help however capable the model, so the lift on covered pairs is the
reading that speaks to the model's ability to use evidence, and the two are
always shown together.

### Fills in as rungs finish

Only rungs with run files are loaded. The 405B rung runs on a cluster; the
tables and figures redraw when its files arrive.

## Run

    BKG_MODE=real BKG_FAMILIES=main BKG_MODELS=llama3.2:1b,llama3.2:3b,llama3.1:8b jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb
    jupyter nbconvert --to notebook --execute grounding_experiment/04_scaling/notebook.ipynb

`capacity_table` and `scaling_table` are in `src/biokgsuite/grounding/analysis.py`;
figures in `src/biokgsuite/viz/grounding.py`.
