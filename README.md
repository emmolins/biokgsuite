# BioKGSuite

BioKGSuite is a framework for evaluating biomedical knowledge graphs (KGs) on 
seven quality dimensions (coverage, annotation accuracy, trustworthiness, 
topology, stability, task performance, and generalisation) through 19 metrics 
computed from each KG's released files and task-specific reference data. This
repository is the code behind the paper *Knowledge graph choice is a design 
decision: multidimensional evaluation of biomedical knowledge graphs with 
BioKGSuite* (Molins, Holmes, Orfanoudaki, Ktena, and González, 2026).

Interactive results: https://emmolins.github.io/biokgsuite/ (built from the result tables by `tools/build_dashboard.py`).

## Layout

| Path | Contents |
|---|---|
| `config/` | `kgs.yaml` (the six graphs, their files and identifier schemes), `params.yaml` (every analysis parameter, by notebook), `kg_snapshots.yaml` (content date per graph), `display_items.yaml` (which output is which manuscript item) |
| `data/` | Inputs. `MANIFEST.md` says how each was obtained; `PROVENANCE.md` (generated) records release, retrieval date and checksum of every file on disk |
| `src/biokgsuite/` | The package the notebooks call: `io`, `graph`, `gold`, `grounding`, `stats`, `viz`, `models` |
| `biokgsuite_dimensions/` | Part I: notebooks `01` to `07`, one per dimension, and `00_summary`, which assembles the dimension profiles (Figs 2 and 3). `CONVENTIONS.md` is the structure they share, `FIGURE_STYLE.md` the figure style |
| `grounding_experiment/` | Part II: five stage notebooks on the post-cutoff approval set (Figs 5 to 8) |
| `supplementary/` | Embedding cross-checks of the heuristic results (needs torch) |
| `results/` | Everything the notebooks write: `figures/`, `tables/`, `checkpoints/`. `MANIFEST.md` (generated) lists every file with its checksum; `TRACEABILITY.md` maps manuscript items to files |
| `tests/` | Unit tests of the metrics, the reference-set builder and the grounding package |
| `tools/` | Fetch, build and record scripts; `tools/README.md` groups them |
| `hpc/` | Slurm scripts for the long runs |
| `paper/figures/` | The display items drawn or assembled outside the notebooks: the Fig 1 and Fig 5 schematics, Fig 4 (`make_fig4.py`), Fig 7 (`make_fig7.py`), Fig 8 and Supplementary Fig S3. Everything else is read from `results/figures/` |

## What is in git and what is fetched

The graphs (25 GB), ontologies and Open Targets datasets are not committed;
`.gitignore` keeps only their release notes, checksums and snapshot records,
and `data/MANIFEST.md` gives the fetch commands. Committed inputs are the
derived bridges, the reference sets and the post-cutoff approval set with its
vendored source, about 60 MB. Per-graph views under
`results/checkpoints/07_kg_views/` (300 MB) are rebuilt by
`tools/build_post_cutoff_approvals.py --extract-views`.

## Environment

Python 3.11. `environment.yml` is the reference set and `requirements.txt`
mirrors it for pip; `requirements-embeddings.txt` adds torch for
`supplementary/`.

```bash
bash setup_local.sh                 # creates .venv from requirements.txt (add --embeddings for torch)
source .venv/bin/activate
python -m pytest tests -q           # no install step: pyproject.toml puts src/ on the path
```

Keys for the hosted models in Part II go in `.env` at the repository root
(copy `.env.example`); the file is gitignored and read by
`grounding_experiment/02_run`.

## Reproducing the results

```
1. python tools/fetch_kgs.py --all; python tools/build_crosswalks.py; bash tools/fetch_ctd.sh
   python tools/provenance.py                                   -> data/, data/PROVENANCE.md
2. biokgsuite_dimensions/01_coverage ... 07_generalization      -> results/, one checkpoint each
3. biokgsuite_dimensions/00_summary                             -> Figs 2 and 3, Tables S3 and S4
4. python tools/build_post_cutoff_approvals.py --extract-views <kgs>; then --kgs all --figure
   python tools/build_distractor_pool.py; python tools/build_entrez_symbol_map.py
5. grounding_experiment/01_benchmark                            -> benchmark and selection audit (Table S5)
6. grounding_experiment/02_run            [API keys, hours]     -> results/tables/G02_run/runs/
7. grounding_experiment/03_analysis, 04_scaling, 05_diagnostic  -> Fig 6, Figs S4 and S5, Table S6
8. supplementary/embeddings               [GPU]                 -> Table S7
9. python tools/results_manifest.py; python tools/regenerate_figs.py
   python paper/figures/make_fig4.py; python paper/figures/make_fig7.py   -> results/MANIFEST.md, Figs 4 and 7
```

Step 4 is committed: `data/gold_standards/post_cutoff_approvals/` rebuilds
byte-identically (its `PROVENANCE.txt` carries the source SHA-256). The
grounding arms' views also carry the grounded flag, so they are re-extracted
after any change to the dossier rule. Only step 6 is expensive, and it is
idempotent per (family, model, seed) file. On a Slurm cluster the whole order
is `hpc/README.md`, one job script per step. Each notebook is run from its own
folder:

```bash
cd biokgsuite_dimensions/05_stability
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=86400 notebook.ipynb
```

The notebooks are committed with the outputs of the run behind the manuscript,
so every number in the tree comes from one run. Figures redraw from the
checkpoints without re-running anything: `python tools/regenerate_figs.py --only 05`.

To trial the notebooks without loading the whole panel (MATRIX alone is 81 M
edges), write the graphs to load into a `.trial_kgs` file at the repository
root, for example `echo primekg,hetionet > .trial_kgs`, or set `BKG_KGS` in
the environment before the kernel starts. Every notebook and tool then sees
only those graphs and prints a trial-run notice, and the checkpoint records
the run as a subset. `.trial_kgs` is gitignored, so a clone on the cluster
always runs the full six.

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
it differs between graphs, is in `data/MANIFEST.md`. Superseded work is kept
at the tag `v0.9-legacy` and is not read by anything on `main`.

## Citation

If you use BioKGSuite, please cite the manuscript above and the software
release (`CITATION.cff`):

```bibtex
@software{biokgsuite2026,
  author  = {Molins, Emily and Holmes, Chris and Orfanoudaki, Agni and Ktena, Ira and Gonz{\'a}lez, Javier},
  title   = {BioKGSuite: a multidimensional evaluation framework for biomedical knowledge graphs},
  version = {1.0.2},
  year    = {2026},
  url     = {https://github.com/emmolins/biokgsuite}
}
```

## Licence

MIT (`LICENSE`). The data sources carry their own licences, listed per file in
`data/PROVENANCE.md`.
