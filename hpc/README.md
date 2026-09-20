# Running everything on the cluster

Target: Oracle Cloud HPC (Ashburn), Slurm, partition `def`. CPU nodes have 96
cores and 125 GB or more; GPU nodes carry H100 80 GB cards. The repo lives at
`~/biokgsuite-bib` on the cluster and every script below assumes that path
(`REPO_DIR` overrides it).

The runs, in the order they depend on each other:

| Job | Script | Needs | Takes | Produces |
|---|---|---|---|---|
| A. Dimensions 01 to 07, then 00 | `run_main.sbatch` | 256 GB CPU node | 6 to 12 h | `results/`, executed notebooks, Fig 1 |
| B. Grounding benchmark | `run_grounding_analysis.sbatch STAGE=benchmark` | CPU | 10 min | `G01_benchmark.pkl`, Fig G01 |
| C. Grounding runs, hosted models | `run_grounding_api.sbatch` x 9 | `.env`, B | 3 to 6 h each | `results/tables/G02_run/runs/` |
| D. Llama ladder, 1B to 70B | `run_llm_ladder.sbatch` x 12 | 1 H100, B | 1 to 8 h each | more run files |
| E. Llama 405B | `run_llm_405b.sbatch` x 3 | 4 H100, B | up to 24 h each | more run files |
| F. Grounding analysis 02 to 05 | `run_grounding_analysis.sbatch` | C (and D, E if wanted) | 30 min | Figs G03 to G05 |
| G. Embedding supplementary | `run_embeddings.sbatch` | 1 H100, HF login | 2 to 4 h | Figs S1 to Sn |
| H. Stability supplementary | `run_stability_embeddings.sbatch` | 1 H100 | hours | supplementary tables |

A and B are independent of each other. C, D and E need B and can all run at
once. F runs when the run files you want are in. G and H are independent of
everything else.

## 0. Cheat sheet

```bash
ssh ashburn-login                       # log in
squeue -u $USER                         # my jobs (ST = PD pending, R running)
tail -f logs/<name>-<jobid>.out         # follow a job
scancel <jobid>                         # kill one job;  scancel -u $USER  kills them all
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS   # after it ends: did it finish, how much memory
srun --partition=def --cpus-per-task=8 --mem=64G --time=4:00:00 --pty bash   # interactive shell on a node
```

## 1. On the laptop: commit and push

The repo has no remote yet. Once, from the laptop (a stale lock from an
interrupted commit has to go first):

```bash
cd ~/Desktop/biokgsuite-bib
rm -f .git/index.lock && find .git/objects -name 'tmp_obj_*' -delete
git add -A
git commit -m "Identifier resolution, reference sets, grounded stage, HPC scripts"
gh repo create emmolins/biokgsuite-bib --private --source=. --remote=origin --push
```

Every later change is the usual pair:

```bash
cd ~/Desktop/biokgsuite-bib && git push origin main          # laptop
cd ~/biokgsuite-bib && git pull --ff-only                    # cluster
```

`.env`, `.trial_kgs`, the graphs under `data/` and the loader caches are
gitignored, so a pull never carries a key, never carries the two-graph trial
switch, and never moves data. That is what you want: the cluster runs the
full panel because `.trial_kgs` is absent there, and every job script refuses
to start if it finds one.

## 2. On the cluster, once: clone, environment, data, keys

```bash
ssh ashburn-login
git clone https://github.com/emmolins/biokgsuite-bib.git ~/biokgsuite-bib
cd ~/biokgsuite-bib
bash hpc/setup.sh            # uv venv at .venv, requirements.txt + the embedding stack; ~10 min
```

This repo has its own `.venv`; do not activate another project's. The job
scripts activate it themselves.

The graphs are not in git. The fast way is to let the cluster fetch them
from the publishers at the pinned releases and check every derived file
against the record committed in `data/<kg>/<KG>_RELEASE.txt`, so the cluster
is known to hold the laptop's bytes:

```bash
sbatch hpc/fetch_data.sbatch            # primekg, openbilink, matrix; 20 to 40 minutes
tail -f logs/fetch-<jobid>.out          # must end with every line OK and DONE
```

Hetionet, BioKG and DRKG are small and come from the laptop with rsync (DRKG
is fetchable too, but its publisher posts no checksum). From the laptop:

```bash
rsync -avz --partial --progress --exclude '_cache_matrix_*' --exclude 'raw' \
      --include='/hetionet/***' --include='/biokg/***' --include='/drkg/***' --include='/crosswalks/***' \
      --include='/gold_standards/***' --include='/*.md' --exclude='*' \
      ~/Desktop/biokgsuite-bib/data/ ashburn-login:biokgsuite-bib/data/
scp ~/Desktop/biokgsuite-bib/.env ashburn-login:biokgsuite-bib/.env
```

Pushing everything from the laptop instead (about 25 GB) also works and is
the same command without the include lines; the MATRIX loader caches are
left behind because their file name is keyed on the laptop's absolute paths
and would not be reused, and `data/raw` holds only the fetch snapshots the
notebooks never read.

Then on the cluster: `du -sh ~/biokgsuite-bib/data/*/` should show matrix
about 12 GB, openbilink 1.7 GB, primekg 0.9 GB, and `cat .env` should show the
three keys. The first job that loads MATRIX rebuilds its cache on the cluster
(20 minutes or so, once).

Only for the embedding supplementary (G): `huggingface-cli login` with your
token, and accept the EmbeddingGemma licence in a browser.

## 3. Submit

Everything from `~/biokgsuite-bib`. Each `sbatch` prints a job id; the log is
`logs/<name>-<jobid>.out`.

```bash
# A. dimensions 01-07 then 00, one CPU job
sbatch hpc/run_main.sbatch

# B. grounding benchmark; wait for it to finish (squeue) before C, D, E
sbatch --export=ALL,STAGE=benchmark hpc/run_grounding_analysis.sbatch

# Rate limits. A run that stops with "incomplete (rate-limited) arms" in the
# log kept every completed arm; resubmit the same (model, seed) and it finishes
# the file. The error column of the runs file says which limit was hit: a
# per-minute limit clears between retries, a per-day quota does not, and the
# client now aborts at once on a per-day quota. If Gemini reports
# GenerateRequestsPerDay..., the key is on a free tier: enable billing on the
# Google AI Studio project (or use a paid key) before resubmitting, or the job
# will stop again after the same number of calls the next day.
#
# C. hosted models: one job per (model, seed). Same --job-name per provider
#    plus --dependency=singleton in the script means one job per provider at
#    a time (its rate limit), three providers in parallel.
for S in 0 1 2; do
  sbatch --job-name=ground-gpt    --export=ALL,MODEL=gpt-4.1-mini,SEED=$S                 hpc/run_grounding_api.sbatch
  sbatch --job-name=ground-gemini --export=ALL,MODEL=gemini:gemini-3.5-flash-lite,SEED=$S hpc/run_grounding_api.sbatch
  sbatch --job-name=ground-groq   --export=ALL,MODEL=groq:llama-3.3-70b-versatile,SEED=$S hpc/run_grounding_api.sbatch
done

# D. Llama ladder on single H100s, one job per (model, seed)
for M in llama3.2:1b llama3.2:3b llama3.1:8b llama3.3:70b; do
  for S in 0 1 2; do sbatch --export=ALL,MODEL=$M,SEED=$S hpc/run_llm_ladder.sbatch; done
done

# E. Llama 405B on 4 H100s, one seed at a time; read the pre-flight notes in the script first
sbatch --export=ALL,SEED=0 hpc/run_llm_405b.sbatch

# F. once every run you want is in: 02 (reloads all files), 03, 04, 05
sbatch hpc/run_grounding_analysis.sbatch

# G, H. supplementaries
sbatch hpc/run_embeddings.sbatch
sbatch hpc/run_stability_embeddings.sbatch
```

Seeds are 0, 1, 2 (`config/params.yaml`, `grounding.seeds`). The model names
are the keys of `grounding.models` and `grounding.capacity_models` there; a
name not in the config is a config edit, not a command-line one.

## 4. What to check while it runs

`run_main`: the log prints one banner per notebook; 01 and 02 are the MATRIX
memory peak, 04 (Louvain) and 06 the slow ones. `sacct ... MaxRSS` after the
job says whether 256 GB was needed.

Grounding jobs: the first lines say `MODE=real models=[...] families=[...]`
and `loading primekg ...`; a `tqdm` line then appears every 30 s. Watch for
`incomplete (rate-limited) arms` at the end: that job stopped early on eight
consecutive rate-limit errors. Resubmit the same command; the runner keeps the
file's complete arms and runs only the missing ones. `pos_has_evidence agrees
with grounded_<kg>` must print on every run job; a `RuntimeError` there means
the dossier and the status builder disagree and nothing should be analysed
until that is understood.

`ls results/tables/G02_run/runs/` shows one `<family>__<model>__seed<k>.csv`
per finished (family, model, seed); C alone gives 18 files, D adds 12 and E
adds 3 (`main` only for the ladder).

## 5. Pull results back to the laptop

When `squeue -u $USER` is empty, from the laptop:

```bash
rsync -avz ashburn-login:biokgsuite-bib/results/ ~/Desktop/biokgsuite-bib/results/
rsync -avz ashburn-login:biokgsuite-bib/biokgsuite_dimensions/ ~/Desktop/biokgsuite-bib/biokgsuite_dimensions/
rsync -avz ashburn-login:biokgsuite-bib/grounding_experiment/ ~/Desktop/biokgsuite-bib/grounding_experiment/
rsync -avz ashburn-login:biokgsuite-bib/supplementary/ ~/Desktop/biokgsuite-bib/supplementary/
```

That brings the figures, tables, checkpoints, run files and the executed
notebooks with their outputs. Then on the laptop: `python tools/results_manifest.py`,
`python tools/regenerate_figs.py`, commit, push. Leave the laptop's
`.trial_kgs` where it is; it only affects the dimension notebooks when they
are run there, and the checkpoints just pulled from the cluster carry
`full_panel: true`.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `module: command not found` | login environment not loaded | `source /etc/profile.d/modules.sh` first |
| `FATAL: ... is not biokgsuite-bib` | `REPO_DIR` points elsewhere | `export REPO_DIR=~/biokgsuite-bib` or run from that directory |
| `FATAL: .trial_kgs is present` | the laptop's trial switch was copied | `rm ~/biokgsuite-bib/.trial_kgs` |
| `FATAL: no .env` | keys not copied | `scp` step in section 2 |
| `FATAL: run grounding_experiment/01_benchmark first` | B has not run | submit B, wait, resubmit |
| grounding job ends with `incomplete (rate-limited) arms` | provider rate limit | resubmit the same command; only the missing arms rerun |
| `Invalid partition` | partition name differs | `sinfo`, then edit `--partition=` |
| `cuda is False` | CUDA module missing on that node | the script tries `module load cuda`; ask the admins for the module name |
| `GatedRepoError` for Gemma | licence not approved yet | wait for the approval mail, or set `GEMMA_ENABLED = False` in the notebook |
| job killed at the wall time | MATRIX or a long ladder rung | resubmit; every notebook and every run file resumes from what finished |

## legacy/

`run_llm_ollama.sbatch` and `run_llm_ollama_76.sbatch` drove the 2026-09 runs
on the superseded 116-pair and 76-pair sets. Kept for the record.
