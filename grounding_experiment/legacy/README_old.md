# Part II · KG grounding experiment

One expensive runner, four cheap readers.

| Stage | Notebook | Produces |
|---|---|---|
| `01_construction/` | 253→116 prospective set, leakage filter | the benchmark |
| `02_bias_audit/` | kept vs dropped, Fisher OR + Mann–Whitney | **Fig 2** |
| `03_run/` | prompts, LLM clients, KG dossiers, ranking loop | `runs/09_big_*.csv` |
| `04_analysis/` | paired bootstrap, coverage-controlled, multiple-comparison correction | **Fig 3** |
| `05_scaling/` | Llama 3.x ladder, memorisation confound | **Fig 4** |
| `06_diagnostic/` | which pairs are hard and why | **Fig 5** |

`03_run/` is the only stage needing GPUs, API keys and wall time. It has a
`MODE = 'mock'` switch to exercise the pipeline without spending calls.
Everything else reads `runs/`.

`lib/` holds the four modules shared across these stages.
