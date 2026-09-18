"""
Same-family (Llama 3.x) capability-scaling sweep for the LLM x KG study.

Question: does the KG lift (MRR_kg - MRR_nokg) grow, shrink, or stay flat as the
model scales from 1B -> 405B parameters?  A within-family ladder isolates the
capability axis (architecture / training recipe held roughly fixed) in a way the
cross-family GPT/Gemini/Llama comparison cannot.

Memorization confound: bigger/newer models memorise more drug-disease pairs,
inflating the no-KG baseline and mechanically shrinking apparent lift.  So every
figure splits the KG arm into *covered* (the graph actually contains the pair --
memorization can't be the whole story) vs *all arms pooled*.

Single source of truth: notebooks import `compute_scaling` + the `fig_*`
functions; nothing is hard-coded.  Robust to having only a subset of the ladder
run so far (fills in as HPC jobs finish).
"""
import os, glob, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
try:
    from src import figstyle
except ImportError:
    import figstyle
figstyle.apply()
PALETTE, KG, GEN_MARK = figstyle.PALETTE, figstyle.KG, figstyle.GEN_MARK

# ---- the Llama 3.x ladder: raw model tag -> (params in B, generation, nice label) ----
# NB the ladder spans three generations (no clean single-gen 5-point ladder exists in Llama):
#   1B/3B = Llama 3.2 · 8B = Llama 3.1 · 70B = Llama 3.3 · 405B = Llama 3.1
# Generation is carried through so figures can mark it and the 3.1 spine (8B->405B) reads cleanly.
LADDER = {
    'llama3.2:1b':   (1.0,   '3.2', 'Llama-3.2-1B'),
    'llama3.2:3b':   (3.0,   '3.2', 'Llama-3.2-3B'),
    'llama3.1:8b':   (8.0,   '3.1', 'Llama-3.1-8B'),
    'llama3.1:70b':  (70.0,  '3.1', 'Llama-3.1-70B'),   # optional clean-spine 70B
    'llama3.3:70b':  (70.0,  '3.3', 'Llama-3.3-70B'),   # the 70B already run
    'llama3.1:405b': (405.0, '3.1', 'Llama-3.1-405B'),
}
GEN_MARK = {'3.1': 'o', '3.2': 's', '3.3': 'D'}
C_ALL = '#9AA0A6'    # pooled (all 3 KG arms)
C_COV = '#1B9E77'    # covered arm only
C_NOKG = '#C44E52'   # no-KG baseline


def load_llama_runs(runs_dir):
    """Combine all 09_big_*.csv, keep parsed rows, restrict to Llama-ladder models."""
    fs = glob.glob(os.path.join(runs_dir, '09_big_*.csv'))
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df = df[(df['parsed'] == 1) & (df['model'].isin(LADDER))].copy()
    df['params'] = df['model'].map(lambda m: LADDER[m][0])
    df['gen'] = df['model'].map(lambda m: LADDER[m][1])
    df['label'] = df['model'].map(lambda m: LADDER[m][2])
    return df


def _covmap(cov_csv):
    cov = pd.read_csv(cov_csv)
    cvd = cov.groupby(['disease_name', 'kg'])['covered'].max().reset_index()
    return {(r.disease_name, r.kg): int(r.covered) for r in cvd.itertuples()}


def _boot_ci(vals, n_boot=2000, seed=0):
    """Percentile CI from a cluster (per-disease) bootstrap of a difference vector."""
    if len(vals) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    means = [rng.choice(vals, len(vals), replace=True).mean() for _ in range(n_boot)]
    return tuple(np.percentile(means, [2.5, 97.5]))


def compute_scaling(df, cov_csv):
    """Per-model scaling table: MRR on each arm + KG lift (pooled and covered-only), with bootstrap CIs.

    Lift is computed per disease (mean KG-arm RR - no-KG RR) then averaged, so the CI is a
    cluster bootstrap over the 116 diseases — the same unit used elsewhere in the study.
    """
    cm = _covmap(cov_csv)
    df = df.copy()
    df['covered'] = df.apply(
        lambda r: cm.get((r.disease, r.kg), np.nan) if r.condition == 'kg' else np.nan, axis=1)
    rows = []
    for m, sub in df.groupby('model'):
        params, gen, label = LADDER[m]
        nokg = sub[sub.condition == 'no_kg'].groupby('disease')['reciprocal_rank'].mean()
        kg_all = sub[sub.condition == 'kg'].groupby('disease')['reciprocal_rank'].mean()
        kg_cov = sub[(sub.condition == 'kg') & (sub.covered == 1)].groupby('disease')['reciprocal_rank'].mean()
        lift_all = (kg_all - nokg).dropna().values
        lift_cov = (kg_cov - nokg).dropna().values
        ci_all, ci_cov = _boot_ci(lift_all), _boot_ci(lift_cov)
        rows.append(dict(
            model=m, label=label, params=params, gen=gen,
            mrr_nokg=nokg.mean(), mrr_kg=kg_all.mean(), mrr_kg_cov=kg_cov.mean(),
            lift=lift_all.mean(), lift_lo=ci_all[0], lift_hi=ci_all[1],
            lift_cov=lift_cov.mean(), lift_cov_lo=ci_cov[0], lift_cov_hi=ci_cov[1],
            n_dis=len(lift_all)))
    out = pd.DataFrame(rows).sort_values('params').reset_index(drop=True)
    return out


# ----------------------------- figures (matplotlib, shared house style) -----------------------------
KG_ORDER = ['primekg', 'drkg', 'biokg']
KG_LABEL = {'primekg': 'PrimeKG', 'drkg': 'DRKG', 'biokg': 'BioKG'}


def _logx(ax, params):
    ax.set_xscale('log')
    ax.set_xticks(sorted(set(params)))
    ax.set_xticklabels([f'{int(p)}B' for p in sorted(set(params))])
    ax.minorticks_off()
    ax.set_xlabel('model size (log scale)')


def fig_lift_vs_size(scaling):
    """Headline: KG lift vs model size — pooled vs covered arm, 95% cluster-bootstrap CI."""
    s = scaling.sort_values('params')
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.axhline(0, color='#B4B2A9', ls='--', lw=0.8, zorder=1)
    for arm, col, mk, lo, hi, val in [
        ('covered arm', PALETTE['covered'], 'D', s.lift_cov_lo, s.lift_cov_hi, s.lift_cov),
        ('pooled (all 3 KGs)', PALETTE['pooled'], 'o', s.lift_lo, s.lift_hi, s.lift),
    ]:
        ax.vlines(s.params, lo, hi, color=col, lw=1, zorder=2)
        ax.plot(s.params, val, marker=mk, ls='-', color=col, lw=1.6,
                mfc='white', mec=col, mew=1.4, ms=8, label=arm, zorder=3)
    _logx(ax, s.params)
    ax.set_ylabel('MRR lift  (KG − no-KG)')
    figstyle.title(ax, 'KG lift vs model size',
                   'Mean reciprocal-rank lift over no-KG  ·  pooled vs covered arm  ·  95% CI')
    ax.legend(loc='lower right', handlelength=1.8)
    fig.tight_layout()
    return fig


def kg_by_model_table(df):
    """Per (model, KG) MRR on the KG arm + the per-model no-KG baseline. Tidy table."""
    d = df.copy()
    present = sorted([m for m in LADDER if m in set(d.model)], key=lambda m: LADDER[m][0])
    rows = []
    for m in present:
        nokg = d[(d.model == m) & (d.condition == 'no_kg')]['reciprocal_rank'].mean()
        rec = {'model': LADDER[m][2], 'params': LADDER[m][0], 'no_KG': nokg}
        for kg in KG_ORDER:
            sub = d[(d.model == m) & (d.condition == 'kg') & (d.kg == kg)]
            rec[KG_LABEL[kg]] = sub['reciprocal_rank'].mean()
        rows.append(rec)
    return pd.DataFrame(rows)


def fig_kg_by_model(df):
    """Grouped bars: each KG's MRR (+ no-KG baseline) across the Llama ladder.
    Shows whether one graph drives the lift and how each KG's benefit scales with model size."""
    t = kg_by_model_table(df).sort_values('params')
    labels = list(t['model'])
    arms = [('no-KG', 'no_KG'), ('PrimeKG', 'PrimeKG'), ('DRKG', 'DRKG'), ('BioKG', 'BioKG')]
    x = np.arange(len(labels))
    w = 0.2
    fig, ax = plt.subplots(figsize=(9.4, 4.8))
    for i, (lab, col) in enumerate(arms):
        ax.bar(x + (i - 1.5) * w, t[col].values, w, label=lab,
               color=KG[lab], edgecolor='white', linewidth=0.6, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('MRR on the 116-set')
    ax.set_ylim(0, None)
    figstyle.title(ax, 'MRR by knowledge graph and model size',
                   'Mean reciprocal rank per KG vs the no-KG baseline, across the Llama ladder')
    ax.legend(ncol=4, loc='upper left', columnspacing=1.3, handlelength=1.1)
    fig.tight_layout()
    return fig


def fig_mrr_vs_size(scaling):
    """Diagnostic: prior-only vs KG-grounded MRR across size — a flat no-KG line rules out memorization."""
    s = scaling.sort_values('params')
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for arm, col, mk, ls, val in [
        ('KG arm (covered)', PALETTE['covered'], 'D', '-', s.mrr_kg_cov),
        ('KG arm (pooled)', PALETTE['pooled'], 'o', '-', s.mrr_kg),
        ('no-KG (prior only)', PALETTE['nokg'], 's', '--', s.mrr_nokg),
    ]:
        ax.plot(s.params, val, marker=mk, ls=ls, color=col, lw=1.6,
                mfc='white', mec=col, mew=1.4, ms=7.5, label=arm, zorder=3)
    _logx(ax, s.params)
    ax.set_ylabel('MRR on the 116-set')
    figstyle.title(ax, 'Prior-only vs KG-grounded MRR by model size',
                   'no-KG, pooled KG and covered-arm MRR across the Llama ladder')
    ax.legend(loc='lower right', handlelength=2.2)
    fig.tight_layout()
    return fig
