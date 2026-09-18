"""
"Which drug-disease pairs rank poorly?" analysis for the LLM x KG study.

compute_difficulty() builds a per-disease difficulty table from the run data
(joined to drug identity, therapeutic area, KG coverage). The three plotting
functions visualise it. Dataset-agnostic given the same run schema.
"""
import os, numpy as np, pandas as pd
import matplotlib                       # don't force a backend (keeps %matplotlib inline working)
import matplotlib.pyplot as plt
from .headline_figures import load_runs, MODELS, MCOL
try:                                    # shared house style (consistent across 09/09a/09b/09c)
    from . import figstyle
except ImportError:
    import figstyle
figstyle.apply()

_MC = [m.split('-')[0] for m in MODELS]          # short keys: GPT, Gemini, Llama


import re as _re

def _modality(active, drug):
    s = (str(active) + ' ' + str(drug)).lower()
    if _re.search(r'cabtagene|leucel|car-t', s): return 'Cell / gene'
    if 'vaccine' in s: return 'Vaccine'
    if _re.search(r'mab\b|umab|zumab|ximab|cept\b', s): return 'Antibody'
    if _re.search(r'tide\b|glutide|parin', s): return 'Peptide'
    return 'Small molecule'

def compute_difficulty(df, gold_tsv, audited_csv=None, pairs_csv=None):
    """Per-disease difficulty on the KG arm, pooled over models + seeds, joined to
    drug identity, therapeutic area, KG coverage, drug age (orig approval year) and modality."""
    kg = df[df.condition == 'kg']
    g = pd.read_csv(gold_tsv, sep='\t')
    meta = g.groupby('disease_name').agg(
        drug=('drug_name', lambda s: ' / '.join(sorted(set(s)))),
        area=('therapeutic_area', 'first'), doid=('disease_id', 'first')).reset_index()
    per = kg.groupby('disease').agg(mrr=('reciprocal_rank', 'mean'), hit1=('hit@1', 'mean'),
                                    conf=('pos_confidence', 'mean'), n=('reciprocal_rank', 'size')).reset_index()
    for m, k in zip(MODELS, _MC):
        per = per.merge(kg[kg.Model == m].groupby('disease')['reciprocal_rank'].mean().rename(k),
                        on='disease', how='left')
    per['frac_consistent'] = kg.assign(c=(kg.pos_agreement == 'consistent').astype(int)) \
        .groupby('disease')['c'].mean().reindex(per['disease']).values
    per = per.merge(meta, left_on='disease', right_on='disease_name', how='left')
    if audited_csv and os.path.exists(audited_csv):
        a = pd.read_csv(audited_csv); a['n3'] = a.eval_primekg + a.eval_drkg + a.eval_biokg
        per['n_kg'] = per['doid'].map(a.groupby('DOID_ID')['n3'].max())
    # drug age + modality from the candidate file (joined on DOID)
    if pairs_csv is None:
        pairs_csv = os.path.join(os.path.dirname(gold_tsv), 'pairs_with_ids_v2.csv')
    if os.path.exists(pairs_csv):
        p = pd.read_csv(pairs_csv)
        p['orig'] = pd.to_numeric(p.get('OrigYear'), errors='coerce')
        p['modality'] = p.apply(lambda r: _modality(r.get('Active'), r.get('Drug')), axis=1)
        p['doidn'] = p['DOID_ID'].map(lambda x: None if pd.isna(x) else str(x).split(':')[-1].lstrip('0'))
        pm = p.dropna(subset=['doidn']).drop_duplicates('doidn').set_index('doidn')
        dn = per['doid'].map(lambda x: None if pd.isna(x) else str(x).split(':')[-1].lstrip('0'))
        per['orig'] = dn.map(pm['orig']); per['modality'] = dn.map(pm['modality']); per['year'] = dn.map(pm['YearSet'])
    per['model_spread'] = per[_MC].max(axis=1) - per[_MC].min(axis=1)
    return per.sort_values('mrr').reset_index(drop=True)


def cross_model_corr(per):
    return per[_MC].corr().values[np.triu_indices(len(_MC), 1)].mean()


def fig_hardest_pairs(per, n=20):
    """Top-n hardest pairs: pooled mean MRR (bar) + per-model MRR (dots)."""
    h = per.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.4, 0.34 * n + 1.4))
    y = np.arange(len(h))
    ax.barh(y, h['mrr'], color=figstyle.PALETTE['grid'], height=0.64, zorder=1)
    for m, k in zip(MODELS, _MC):
        ax.scatter(h[k], y, s=26, color=MCOL[m], zorder=3, edgecolor='white', lw=0.4, label=m)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{r.disease[:40]}\n({str(r.drug)[:26]} · {r.area})' for r in h.itertuples()], fontsize=6.8)
    ax.set_xlim(0, 1.0); ax.set_xlabel('Mean MRR on KG arm   (grey = pooled · dots = per model)')
    ax.legend(loc='lower right', fontsize=8, title='per-model MRR', title_fontsize=7.5)
    figstyle.title(ax, f'The {n} hardest-to-rank drug–disease pairs', pad=10)
    fig.tight_layout(); return fig


def fig_difficulty_drivers(per):
    """The model knows: per-disease MRR vs self-reported confidence, coloured by
    share of consistent KG evidence."""
    fig, ax = plt.subplots(figsize=figstyle.FIGSIZE)
    sc = ax.scatter(per['conf'], per['mrr'], c=per['frac_consistent'], cmap=figstyle.SEQ_CMAP,
                    s=46, edgecolor='white', lw=0.4, vmin=0, vmax=1)
    r = per[['conf', 'mrr']].corr().iloc[0, 1]
    ax.text(0.04, 0.96, f'Pearson r = {r:.2f}', transform=ax.transAxes, va='top', fontsize=10,
            color=figstyle.PALETTE['ink'], fontweight='bold')
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02); cb.set_label('share of consistent KG evidence', fontsize=9)
    ax.set_xlabel("model's self-reported confidence (1–5)")
    ax.set_ylabel('mean MRR on KG arm')
    figstyle.title(ax, 'Per-pair MRR vs model confidence', pad=10)
    fig.tight_layout(); return fig


_COVCOL = figstyle.COV_RAMP                              # ordinal teal ramp: 1/2/3 KGs covering

def fig_difficulty_by_area(per, min_n=4):
    """Per-area difficulty: each dot a disease, a clear mean marker per area, areas
    sorted easiest→hardest (top→bottom). One colour — kept deliberately simple."""
    keep = per['area'].value_counts(); keep = keep[keep >= min_n].index
    d = per[per['area'].isin(keep)].copy()
    stat = d.groupby('area')['mrr'].agg(['mean', 'size']).sort_values('mean', ascending=False)
    order = stat.index.tolist()                              # easiest first -> top
    rng = np.random.default_rng(7)
    fig, ax = plt.subplots(figsize=(8.0, 0.62 * len(order) + 1.3))
    ax.axvline(per['mrr'].mean(), color='#c8c8c8', ls='--', lw=1, zorder=1)
    for i, a in enumerate(order):
        y = len(order) - 1 - i; sub = d[d.area == a]
        ax.scatter(sub['mrr'], y + rng.uniform(-0.15, 0.15, len(sub)), s=34, color='#BFD8D2',
                   alpha=0.9, edgecolor='white', lw=0.3, zorder=2)
        mu = sub['mrr'].mean()
        ax.scatter(mu, y, s=150, color='#1F6F5E', marker='D', zorder=4, edgecolor='white', lw=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f'{a}  (n={int(stat.loc[a, "size"])})' for a in order[::-1]], fontsize=9.5)
    ax.set_xlim(0, 1.02); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel('MRR on KG arm   (small dot = disease · large diamond = area mean)')
    ax.text(per['mrr'].mean(), len(order) - 0.5, 'overall mean', fontsize=8, color=figstyle.PALETTE['muted'], ha='center')
    figstyle.title(ax, 'Ranking difficulty by therapeutic area', pad=10)
    fig.tight_layout(); return fig


def fig_intrinsic_difficulty(per):
    """Hard pairs are hard for every model: pairwise per-disease MRR scatters with the
    y=x line and Pearson r. Points hug the diagonal -> difficulty is intrinsic."""
    from scipy import stats as st
    pairs = [(_MC[0], _MC[1]), (_MC[0], _MC[2]), (_MC[1], _MC[2])]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2))
    for ax, (a, b) in zip(axes, pairs):
        d = per[[a, b]].dropna()
        ax.plot([0, 1], [0, 1], color='#bbb', ls='--', lw=1, zorder=1)
        ax.scatter(d[a], d[b], s=34, color=figstyle.PALETTE['covered'], alpha=0.7, edgecolor='white', lw=0.3, zorder=2)
        r = st.pearsonr(d[a], d[b])[0]
        ax.text(0.04, 0.93, f'r = {r:.2f}', transform=ax.transAxes, fontsize=11, fontweight='bold', color=figstyle.PALETTE['ink'])
        ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02); ax.set_aspect('equal')
        ax.set_xlabel(f'{a} — MRR'); ax.set_ylabel(f'{b} — MRR')
    ravg = cross_model_corr(per)
    figstyle.suptitle(fig, f'Per-pair MRR agreement between models   (mean pairwise r = {ravg:.2f})', y=1.02)
    fig.tight_layout(); return fig


def fig_coverage_fair(df, cov_csv):
    """Is the coverage→difficulty effect real or an averaging artifact? Compares, by the
    pair's total KG coverage, the POOLED per-disease MRR (all 3 KG arms, incl. empty ones)
    with MRR on only the arm(s) that actually cover the pair. Flat covered-arm line =>
    the pooled slope is mostly pairs being absent from some graphs."""
    kg = df[df.condition == 'kg'].copy()
    cov = pd.read_csv(cov_csv)
    # per-(disease, kg) covered (max over multiple drugs of a disease), then n_kg per disease
    cvd = cov.groupby(['disease_name', 'kg'])['covered'].max().reset_index()
    covmap = {(r.disease_name, r.kg): int(r.covered) for r in cvd.itertuples()}
    n_kg = cvd.groupby('disease_name')['covered'].sum()
    kg['covered'] = kg.apply(lambda r: covmap.get((r.disease, r.kg), np.nan), axis=1)
    kg['n_kg'] = kg['disease'].map(n_kg)
    kg = kg[kg['n_kg'].isin([1, 2, 3])]
    pooled = kg.groupby('n_kg')['reciprocal_rank'].mean()
    covered = kg[kg.covered == 1].groupby('n_kg')['reciprocal_rank'].mean()
    nokg = df[df.condition == 'no_kg']['reciprocal_rank'].mean()
    x = [1, 2, 3]
    _pool, _cov = figstyle.PALETTE['pooled'], figstyle.PALETTE['covered']
    fig, ax = plt.subplots(figsize=figstyle.FIGSIZE)
    ax.axhline(nokg, color='#bbb', ls=':', lw=1.2); ax.text(3.05, nokg, ' no-KG', va='center', fontsize=8.5, color=figstyle.PALETTE['muted'])
    ax.plot(x, [pooled.get(k, np.nan) for k in x], '-o', color=_pool, lw=1.8, ms=9,
            mfc='white', mec=_pool, mew=1.4, label='pooled (all 3 KG arms)')
    ax.plot(x, [covered.get(k, np.nan) for k in x], '-D', color=_cov, lw=1.8, ms=9,
            mfc='white', mec=_cov, mew=1.4, label='covered arm only')
    for k in x:
        if k in covered: ax.text(k, covered[k] + 0.03, f'{covered[k]:.2f}', ha='center', fontsize=9, color=_cov, fontweight='bold')
        if k in pooled: ax.text(k, pooled[k] - 0.045, f'{pooled[k]:.2f}', ha='center', fontsize=9, color=figstyle.PALETTE['muted'], fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels([f'{k} KG' + ('s' if k > 1 else '') for k in x])
    ax.set_xlim(0.7, 3.5); ax.set_ylim(0, 0.9)
    ax.set_xlabel('number of KGs covering the pair'); ax.set_ylabel('MRR on KG arm')
    ax.legend(loc='lower right')
    figstyle.title(ax, 'MRR by number of KGs covering the pair',
                   "Pooled (all 3 arms) vs covered-arm only — the gap is the empty-arm effect")
    fig.tight_layout(); return fig


def fig_age_effect(per):
    """Drug age does not predict difficulty: original approval year vs MRR."""
    d = per.dropna(subset=['orig', 'mrr'])
    from scipy import stats as st
    r, p = st.pearsonr(d['orig'], d['mrr'])
    s, i = np.polyfit(d['orig'], d['mrr'], 1)
    fig, ax = plt.subplots(figsize=figstyle.FIGSIZE)
    ax.scatter(d['orig'], d['mrr'], s=44, color=figstyle.PALETTE['pooled'], alpha=0.75, edgecolor='white', lw=0.4, zorder=2)
    xl = np.array([d['orig'].min(), d['orig'].max()])
    ax.plot(xl, s * xl + i, color=figstyle.PALETTE['nokg'], lw=1.8, ls='--', zorder=3)
    ax.text(0.04, 0.06, f'Pearson r = {r:+.2f}   (p = {p:.2f}, n.s.)', transform=ax.transAxes,
            fontsize=10.5, color=figstyle.PALETTE['nokg'], fontweight='bold')
    ax.set_xlabel("drug's original approval year   (older ←      → newer)")
    ax.set_ylabel('mean MRR on KG arm'); ax.set_ylim(0, 1.02)
    figstyle.title(ax, 'MRR vs drug approval year', pad=10)
    fig.tight_layout(); return fig


def fig_coverage_effect(per):
    """The axis that actually predicts difficulty: number of KGs covering the pair."""
    d = per.dropna(subset=['n_kg']).copy(); d['n_kg'] = d['n_kg'].astype(int)
    groups = [1, 2, 3]; rng = np.random.default_rng(11)
    fig, ax = plt.subplots(figsize=figstyle.FIGSIZE)
    for i, k in enumerate(groups):
        v = d[d.n_kg == k]['mrr']
        ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v, s=44, color=_COVCOL[k], alpha=0.85,
                   edgecolor='white', lw=0.4, zorder=3)
        ax.plot([i - 0.28, i + 0.28], [v.mean()] * 2, color=figstyle.PALETTE['ink'], lw=2.6, zorder=4, solid_capstyle='round')
        ax.text(i, v.mean() + 0.03, f'{v.mean():.2f}', ha='center', fontsize=9.5, fontweight='bold', color=figstyle.PALETTE['ink'])
    ax.set_xticks(range(3)); ax.set_xticklabels([f'{k} KG' + ('s' if k > 1 else '') + f'\n(n={int((d.n_kg==k).sum())})' for k in groups])
    ax.set_ylim(0, 1.02); ax.set_ylabel('MRR on KG arm  (dot = disease)')
    ax.set_xlabel('number of KGs covering the pair')
    figstyle.title(ax, 'MRR by KG coverage count', pad=10)
    fig.tight_layout(); return fig
