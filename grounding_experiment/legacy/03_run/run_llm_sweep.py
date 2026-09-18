#!/usr/bin/env python
"""Standalone BioKGSuite LLM ranking sweep — one model -> one CSV.

Executes the setup cells + the ranking-loop cell of 09_llm_integration.ipynb.
Cells are located by CONTENT (the run cell is the one defining `run_arm`), not by
hardcoded index, so this survives the notebook being re-ordered. It also forces
real-run mode (LOAD_COMPILED=False, FORCE_RERUN=True) so the loop actually calls
the model instead of loading previously-compiled CSVs.

Usage:
    python eval_notebooks/run_llm_sweep.py --model llama3.2:1b --gold gold_standard_v4_bigset.tsv --seeds 1 --shuffles 2 --out 09_big_llama3_2_1b_s1.csv
"""
import os, sys, json, argparse

ap = argparse.ArgumentParser()
ap.add_argument('--model', required=True, help="e.g. llama3.2:1b | gpt-4.1-mini | gemini:gemini-3.5-flash")
ap.add_argument('--out', help="output CSV name (default 09_ranking_<model>.csv)")
ap.add_argument('--seeds', default=None, help="comma seeds, e.g. 0 or 0,1")
ap.add_argument('--shuffles', default=None)
ap.add_argument('--n', default=None, help='limit N_DISEASES (for a quick smoke test)')
ap.add_argument('--gold', default=None, help='gold standard filename under data/gold_standards/legacy/')
a = ap.parse_args()

os.environ['BKG_MODELS'] = a.model
if a.seeds:    os.environ['BKG_SEEDS'] = a.seeds
if a.shuffles: os.environ['BKG_SHUFFLES'] = a.shuffles
if a.n:        os.environ['BKG_NDISEASES'] = a.n
if a.gold:     os.environ['BKG_GOLD'] = a.gold
safe = a.model.replace(':', '_').replace('/', '_')
os.environ['BKG_OUT_CSV'] = a.out or f'09_ranking_{safe}.csv'

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)                       # so the setup cell resolves repo root
nb = json.load(open('09_llm_integration.ipynb'))
cells = nb['cells']

# locate the ranking-loop cell by content (defines run_arm) — robust to re-ordering
run_idx = next((i for i, c in enumerate(cells)
                if c['cell_type'] == 'code' and 'def run_arm' in ''.join(c['source'])), None)
if run_idx is None:
    sys.exit("FATAL: no cell defines `run_arm` in 09_llm_integration.ipynb — "
             "this is the wrong 09 version (need the MRR ranking-sweep notebook, not the AUROC one).")
setup_idx = [i for i, c in enumerate(cells) if c['cell_type'] == 'code' and i < run_idx]

g = {}
for i in setup_idx:
    exec(''.join(cells[i]['source']), g)

# force real-run mode regardless of the notebook's saved defaults
g['LOAD_COMPILED'] = False
g['FORCE_RERUN'] = True

print(f"[run] model={g.get('MODELS')}  KGS={g.get('KGS')}  gold={getattr(g.get('GOLD_FILE'), 'name', '?')}  "
      f"N={g.get('N_DISEASES')}  seeds={g.get('SEEDS')}  shuffles={g.get('SHUFFLES')}  "
      f"run_cell={run_idx}  -> {os.environ['BKG_OUT_CSV']}")
exec(''.join(cells[run_idx]['source']), g)
print(f"[done] wrote {os.environ['BKG_OUT_CSV']}")
