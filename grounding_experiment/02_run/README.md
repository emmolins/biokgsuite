# G02 · Run

**Computes.** The rank each model gives the true post-cutoff drug in a pool of
eight candidates, for every query, in a `none` arm and one arm per graph, over
three seeds and three shuffles. One row per call, written per (family, model,
seed) file.

**Inputs.** `config/kgs.yaml`, `config/params.yaml` (`grounding`),
`results/checkpoints/G01_benchmark.pkl`, the KG files named in `kgs.yaml`
(real mode only), `data/gold_standards/bridges/` (including `entrez_to_symbol.tsv`,
without which DRKG gene nodes would render as numbers)
**Outputs.** `results/tables/G02_run/runs/<family>__<model>__seed<k>.csv`,
`results/tables/G02_run/worked_examples/`, `results/tables/G02_run/run_summary.csv`,
`results/checkpoints/G02_run.pkl`
**Display item.** None (the worked examples go in the supplement)

## What is measured

| | measured | reported, not scored |
|---|---|---|
| | rank of the positive, reciprocal rank, hit@1/3/5, parsed, error | the model's self-reported `basis`, `evidence_agreement` and `confidence`; whether the positive and how many candidates carried evidence |

A failed parse is a miss (`rank = pool_size + 1`), never dropped.

## Design notes

### The only thing that changes between arms is the evidence

Same pairs, same distractors, same seeds and shuffles in every arm. In a KG
arm each candidate carries that graph's dossier (its gene neighbours with the
graph's own interaction verb, its pathways) and the disease carries a profile
(associated genes, pathways, related diseases). The drug-disease bridge is
excluded structurally: a disease neighbour of a drug is never rendered, and
the target disease is masked in every dossier, so no arm can leak the answer
by construction rather than by a per-graph relation list.

### Candidates are anonymised

Candidates appear as `Drug-1` to `Drug-8`, numbered in the shuffled order, so
the model cannot recall the approval from the drug's name and the `none` arm
measures pharmacological reasoning over the disease alone. Gene names are
shown (`anonymize_genes: false`) because the mechanism is what the evidence
is for.

### Repeats and seeds

Three seeds redraw the distractors; three shuffles reorder the candidates.
Each (seed, pair, shuffle, arm) has its own random stream, so a candidate
order depends on nothing that ran before it and a partial rerun reproduces
the same prompts. The model call seed is derived the same way.

### Idempotent and scoped

A (family, model, seed) file that exists is loaded, not recomputed, unless
`BKG_FORCE=1`. `BKG_MODELS`, `BKG_SEEDS` and `BKG_FAMILIES` scope one launch
to one job. Eight consecutive rate-limit errors abandon the arm for that model
and the checkpoint records it, so the next launch can complete it rather than
silently scoring the failures as misses.

### Memorisation is not controlled here

The approvals postdate the graphs, not the models. The 2023 to 2026 approval
dates are in the benchmark so G04 and the manuscript can stratify by date
against each model's training cutoff; anonymising the candidates limits, but
does not remove, recall.

## Run

    BKG_MODE=mock jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb
    BKG_MODE=real BKG_MODELS=gpt-4.1-mini BKG_SEEDS=0 jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb

Real mode loads the three graphs once (several GB) and makes
`queries x arms x shuffles` calls per model and seed (128 x 4 x 3 = 1,536
for `main`). Keys: `OPENAI_API_KEY`, `GEMINI_API_KEY`, `GROQ_API_KEY`, or
Ollama for local tags.

The loop is `src/biokgsuite/grounding/runner.py`; prompt and parsing
`prompt.py`; model dispatch `clients.py`; dossiers `dossiers.py`.
