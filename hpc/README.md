> **Paths updated 2026-09-16.** `REPO_DIR` now defaults to `$HOME/biokgsuite-bib`
> and notebooks are executed as `biokgsuite_dimensions/<nb>/notebook.ipynb`.
> The old default pointed at a different checkout with an `eval_notebooks/`
> layout, so a submit ran superseded code silently. Both job scripts now
> refuse to start if `REPO_DIR` is not this repo.
>
> `run_all.sbatch` targets `supplementary/embeddings/notebook.ipynb` (repointed 2026-09-18).

# Running nb08 on HPC, the 5-minute playbook

Target cluster: Oracle Cloud Infrastructure HPC (Ashburn).
Scheduler: Slurm. Hardware per compute node: AMD EPYC 9J14 (96 cores), 125 GB RAM, 8× NVIDIA H100 80GB.

## What you'll do

1. **On your Mac**: kill the current Mac kernel + push the repo to GitHub
2. **On HPC**: clone repo, install deps, push KG data over from your Mac
3. **On HPC**: submit one Slurm job
4. **Wait ~2-4 hours**
5. **rsync results back to your Mac**

## Step 1, On your Mac (2 min)

```bash
# In Jupyter: Kernel → Shutdown
# (The Mac was about to spend 6-8 hours encoding Gemma, we don't need it anymore.)

# Push the repo so HPC can clone it
cd ~/biokgsuite
git add . && git commit -m "HPC scripts" && git push
```

If your repo isn't on GitHub yet, skip the push and use rsync in Step 2 instead:
```bash
rsync -avz --exclude .venv --exclude results/cache --exclude '*.ipynb_checkpoints' \
    ~/biokgsuite/ ashburn-login:biokgsuite/
```

## Step 2, On HPC (~10 min, mostly waiting for `uv pip install`)

SSH in and run the setup script:

```bash
ssh ashburn-login   # or however you ssh in
cd ~

# Get the repo
git clone https://github.com/<you>/biokgsuite_dashboard.git biokgsuite
# (or skip if you used rsync in step 1)
cd biokgsuite

# One-time env setup, installs deps via uv
bash scripts/hpc/setup.sh
```

When that finishes, do the **three manual steps** the setup script reminds you about:

```bash
# (a) HuggingFace auth
huggingface-cli login        # paste your hf_... token, say "y" to git credential

# (b) Accept the Gemma license at:
#     https://huggingface.co/google/embeddinggemma-300m
#     (one click in your browser; wait for Google approval email if not already approved)

# (c) Get the KG data, rsync from your Mac:
#     Run this command FROM YOUR MAC (in a separate terminal), NOT on HPC:
```

On your **Mac** (separate terminal):

```bash
rsync -avz --progress ~/biokgsuite/data/ ashburn-login:biokgsuite/data/
```

Expected size: ~30 GB (mostly MATRIX). Transfer time depends on your home internet, typically 15-60 min on a decent connection. The login node has fast network to the compute nodes, so once it's there, the job runs fast.

Quick sanity check on HPC after rsync finishes:

```bash
du -sh ~/biokgsuite/data/*/
# You should see: matrix ~30G, drkg ~400M, primekg ~1G, etc.
```

## Step 3, Submit the Slurm job (~5 sec)

```bash
cd ~/biokgsuite
sbatch scripts/hpc/run_all.sbatch
```

Slurm prints a job ID like `Submitted batch job 1234567`. That's all you need to do, the rest is automated.

## Step 4, Wait, and monitor (~2-4 hours)

```bash
# See your job's status
squeue -u $USER

# Live log
tail -f logs/nb08-1234567.out

# What the log will show, in order:
#   - module load + venv activation
#   - nvidia-smi shows the H100
#   - existing TransE/RotatE single-run training (~20-40 min total across 6 KGs)
#   - Gemma encoding on GPU (~30-60 min for all 6 KGs)
#   - Multi-rerun resampling loop:
#       Rerun 1: TransE+RotatE on all 6 KGs (~30 min)
#       Rerun 2: same, repeated N_RERUNS times
#   - Gemma re-scoring across reruns (~30 min total)
#   - Aggregation + resampled figures
```

If the job dies for any reason (preemption, transient error), just `sbatch` again, the notebook is resumable via per-rerun caches.

## Step 5, Pull results back to your Mac

When `squeue -u $USER` shows your job is gone (= completed):

On your **Mac**:

```bash
rsync -avz ashburn-login:biokgsuite/results/ ~/biokgsuite/results/
rsync -avz ashburn-login:biokgsuite/supplementary/embeddings/notebook.ipynb \
       ~/biokgsuite-bib/supplementary/embeddings/
```

That brings down:
- `results/tables/08_embedding_comparison_resampled.csv`, the long-form numbers (one row per rerun)
- `results/figures/08_embedding_validation/08_resampled_auroc.{pdf,png}` and `08_lift_over_gemma.{pdf,png}`, the resampling figures
- `supplementary/embeddings/notebook.ipynb`, the executed notebook with all cell outputs

Open the notebook locally to see all the inline figures.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `module: command not found` | login env not loaded | `source /etc/profile.d/modules.sh` first |
| `cuda is False` in log | CUDA module not loaded | The sbatch script already tries `module load cuda`, if it fails, ask cluster admin which CUDA module name to use |
| Slurm rejects job: `Invalid partition` | Partition name differs | Check `sinfo` and edit `--partition=` in run_all.sbatch |
| `GatedRepoError` for Gemma | License not yet approved | Either wait for Google's email, or temporarily set `GEMMA_ENABLED = False` in the nb08 Gemma config cell and resubmit |
| Job hits the 24 hr walltime | MATRIX training was slow | Re-`sbatch`, it resumes from cached per-rerun results, only re-doing what's missing |
| `Permission denied` writing to results/ | Repo on read-only filesystem | Don't put repo on a snapshot mount; `~` is fine |
| Want to skip MATRIX (long pole) | n/a | In the nb08 config cell, change `KG_NAMES = list(config['knowledge_graphs'].keys())` to exclude `'matrix'` |

## Costs

The single sbatch job uses:
- 1 H100 GPU x ~2 hr (for Gemma encoding, only if re-encoding)
- 32 CPU cores x ~2-3 hr (for TransE/RotatE training x N_RERUNS reruns)
- 0 dollars (your account is presumably pre-allocated)

After the job completes, the key artifacts are:
- `results/figures/08_embedding_validation/08_resampled_auroc.png`, TransE/RotatE/Gemma per KG with rerun-empirical 95% CIs
- `results/figures/08_embedding_validation/08_lift_over_gemma.png`, per-rerun lift of each trained model over the Gemma name prior

## legacy/

`run_llm_ollama.sbatch` and `run_llm_ollama_76.sbatch` drove the 2026-09 runs on the
superseded 116-pair and 76-pair sets through a runner that no longer exists. Kept for the
record; `run_llm_ladder.sbatch` and `run_llm_405b.sbatch` now launch
`grounding_experiment/02_run/notebook.ipynb` on the post-cutoff approval set.
