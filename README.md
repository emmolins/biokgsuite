# BioKGSuite

A systematic evaluation of six biomedical knowledge graphs for drug
repurposing: seven intrinsic quality dimensions scored on a shared set of
pinned inputs (Part I), and an experiment asking whether each graph's evidence
helps a language model rank drugs that were approved for a new indication after
the graph was built (Part II). Every number and display item in the manuscript
reproduces from this tree plus the fetched, version-stamped inputs.

## Layout

| Path | Contents |
|---|---|
| `config/` | `kgs.yaml` (the six graphs, their files and id schemes), `params.yaml` (every analysis parameter, by notebook), `kg_snapshots.yaml` (content date per graph) |
| `data/` | Inputs. `MANIFEST.md` says how each was obtained, `PROVENANCE.md` (generated) records release, retrieval date and checksum of every file on disk |
| `src/biokgsuite/` | The package the notebooks call: `io`, `graph`, `gold`, `grounding`, `stats`, `viz`, `models` |
| `biokgsuite_dimensions/` | Part I: notebooks `01` to `07`, one per dimension, and `00_summary` (Fig 1). `CONVENTIONS.md` is the skeleton they all follow, `FIGURE_STYLE.md` the house style every figure is drawn in |
| `grounding_experiment/` | Part II: five stage notebooks on the post-cutoff approval set (Figs 2 to 5) |
| `supplementary/` | Embedding cross-checks of the heuristic results (needs torch) |
| `results/` | Everything the notebooks write: `figures/`, `tables/`, `checkpoints/`; `MANIFEST.md` (generated) lists it |
| `tests/` | Unit tests of the metrics, the set builder and the grounding package |
| `tools/` | Fetch, build and record scripts; `tools/README.md` groups them |
| `hpc/` | Slurm scripts for the long runs |
| `paper/` | Display items assembled by `tools/regenerate_figs.py`, numbered as in the manuscript |

## What is in git and what is fetched

The graphs (25 GB), ontologies and Open Targets datasets are not committed;
`.gitignore` keeps only their release notes, checksums and snapshot records.
`data/MANIFEST.md` gives the fetch commands. Committed inputs are the derived
bridges, the reference sets and the post-cutoff approval set with its vendored
source, about 60 MB. Per-graph views under `results/checkpoints/07_kg_views/`
(300 MB) are rebuilt by `tools/build_post_cutoff_approvals.py --extract-views`.

## Environment

Python 3.11. `environment.yml` is the reference set; `requirements.txt`
mirrors it for pip.

```bash
bash setup_local.sh                 # .venv from requirements.txt (add --embeddings for torch)
source .venv/bin/activate
python -m pytest tests -q           # no install step needed: pyproject.toml puts src/ on the path
```

Hosted-model keys for Part II go in `.env` at the repo root (copy
`.env.example`); it is gitignored and read by `grounding_experiment/02_run`.

## Reproduction order

```
1. python tools/fetch_kgs.py --all; python tools/build_crosswalks.py; bash tools/fetch_ctd.sh
   python tools/provenance.py                                   → data/, data/PROVENANCE.md
2. biokgsuite_dimensions/01_coverage ... 07_generalization      → results/, one checkpoint each
3. biokgsuite_dimensions/00_summary                             → Fig 1
4. python tools/build_post_cutoff_approvals.py --extract-views <kgs>; then --kgs all --figure
   (the grounding arms' views also carry the grounded flag; re-extract them after a dossier-rule change)
   python tools/build_distractor_pool.py; python tools/build_entrez_symbol_map.py
5. grounding_experiment/01_benchmark                            → benchmark, selection audit (Fig 2)
6. grounding_experiment/02_run            [API keys, hours]     → results/tables/G02_run/runs/
7. grounding_experiment/03_analysis, 04_scaling, 05_diagnostic  → Figs 3 to 5
8. supplementary/embeddings               [GPU]                 → Figs S1 to Sn
9. python tools/results_manifest.py; python tools/regenerate_figs.py   → results/MANIFEST.md, paper/
```

Step 4 is committed: `data/gold_standards/post_cutoff_approvals/` rebuilds
byte-identically (its `PROVENANCE.txt` carries the source SHA-256). Only step 6
is expensive, and it is idempotent per (family, model, seed) file. On a Slurm
cluster the whole order is `hpc/README.md`, one job script per step. Each
notebook is run from its own folder:

```bash
cd biokgsuite_dimensions/05_stability
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=86400 notebook.ipynb
```

The notebooks are committed with their outputs cleared; the outputs of the
run behind the manuscript are committed with the final results, so that every
number in the tree comes from one run. Figures redraw from the checkpoints
without re-running anything: `python tools/regenerate_figs.py --only 05`.

To trial the notebooks on a laptop without loading the whole panel (MATRIX
alone is 81 M edges), write the graphs to load into a `.trial_kgs` file at the
repo root, for example `echo primekg,hetionet > .trial_kgs`, or set `BKG_KGS`
in the environment before the kernel starts. Every notebook and tool then sees
only those graphs and prints a trial-run notice; the checkpoint records the run
as a subset, so it cannot be mistaken for the full panel. `.trial_kgs` is
gitignored, so a clone on the cluster always runs the full six.

## Knowledge graphs

| KG | Release | Last checked against source |
|---|---|---|
| hetionet | v1.0 | 2026-09-15 |
| primekg | V2.1 (Dataverse `doi:10.7910/DVN/IXA7BM`) | 2026-09-18 |
| drkg | unversioned, 2022-04-19 build | 2026-09-18 |
| openbilink | OpenBioLink2020 v1.0.1 (Zenodo `doi:10.5281/zenodo.5361324`) | 2026-09-18 |
| biokg | v1.0.0 | 2026-09-14 |
| matrix | Every Cure KG, HF `main` @ `c9c60b2c` / `fcd05859` | 2026-09-18 |

`data/PROVENANCE.md` is the authoritative record and is regenerated from the
files on disk; this table is a summary. What backs each version claim, and why
it differs between graphs, is in `data/MANIFEST.md`.

## Housekeeping

`_to_delete/` holds files retired during the 2026-09 clean-up (duplicates,
empty directories, OS files); `bash tools/cleanup.sh --yes` removes it together
with caches. Superseded work is kept, not deleted, under `legacy/` folders
(`data/gold_standards/legacy/`, `grounding_experiment/legacy/`, `tools/legacy/`,
`hpc/legacy/`), each with a README saying what it was and why it was replaced.

## Citation

See `CITATION.cff`. Licence: MIT (`LICENSE`); the data sources carry their own
licences, listed per file in `data/PROVENANCE.md`.
