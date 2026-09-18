"""
Headline LLM x KG figures for notebook 09 (reproducible).

Each function takes the combined run dataframe (and, where needed, the benchmark
summary) and returns a matplotlib Figure. Styling — model colours, KG colours,
footer caption — is shared so the set is visually consistent.

Run/refresh all of them with `regenerate(RUNS_DIR, SUMMARY_CSV, OUT_DIR)`.
"""
import os, glob, numpy as np, pandas as pd
from scipy import stats as st
import matplotlib                      # NB: do not force a backend here — that would
import matplotlib.pyplot as plt        # override %matplotlib inline and hide notebook figures.
from matplotlib import font_manager
# headless scripts (no display) fall back to Agg automatically, so savefig still works.

# ---- shared style ----
_pref = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']
_avail = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in _pref if f in _avail), 'DejaVu Sans')
plt.rcParams.update({'font.family': FONT, 'axes.linewidth': 0.8, 'savefig.dpi': 300, 'figure.dpi': 150,
    'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
try:                                    # shared house style (consistent across 09/09a/09b/09c)
    from . import figstyle
except ImportError:
    import figstyle
figstyle.apply()

MODELS = ['GPT-4.1-mini', 'Gemini-3.1-flash-lite', 'Llama-3.3-70B']
_MSHORT = {'GPT-4.1-mini': 'GPT', 'Gemini-3.1-flash-lite': 'Gemini', 'Llama-3.3-70B': 'Llama'}
MCOL = {m: figstyle.MODEL[_MSHORT[m]] for m in MODELS}     # shared model palette
VCOL = MCOL
VMK = {'GPT-4.1-mini': 'o', 'Gemini-3.1-flash-lite': 's', 'Llama-3.3-70B': 'D'}
KGS = ['primekg', 'drkg', 'biokg']
KLAB = {'primekg': 'PrimeKG', 'drkg': 'DRKG', 'biokg': 'BioKG'}
KCOL = {k: figstyle.KG[KLAB[k]] for k in KGS}              # shared KG palette
RANK_BINS = ['Rank 1', 'Rank 2', 'Rank 3', 'Rank 4–8', 'Missed']
RANK_COL = figstyle.RANK_RAMP                              # shared best->worst ramp
FOOT = '116-set · 3 seeds × 2 shuffles · Open Targets distractors · GPT-4.1-mini + Gemini-3.1-flash-lite + Llama-3.3-70B'
_RAW2NICE = {'gpt-4.1-mini': 'GPT-4.1-mini', 'gemini:gemini-3.1-flash-lite': 'Gemini-3.1-flash-lite',
             'llama3.3:70b': 'Llama-3.3-70B'}


def _footer(fig, text=FOOT):
    fig.text(0.5, 0.005, text, ha='center', fontsize=7.5, color='#9a9a9a')


def load_runs(runs_dir):
    """Combine all 09_big_*.csv, keep parsed rows, normalise model + arm."""
    df = pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(runs_dir, '09_big_*.csv'))], ignore_index=True)
    df = df[df['parsed'] == 1].copy()
    df['Model'] = df['model'].map(_RAW2NICE).fillna(df['model'])
    df['arm'] = np.where(df['condition'] == 'no_kg', 'no_kg', df['kg'])
    return df


def _pdz(df, model, arm):
    return df[(df.Model == model) & (df.arm == arm)].groupby('disease')['reciprocal_rank'].mean()

def _diffs(df, model, kg):
    j = pd.concat([_pdz(df, model, kg).rename('k'), _pdz(df, model, 'no_kg').rename('n')], axis=1).dropna()
    return (j['k'] - j['n']).values

def _lift_ci(df, model, kg):
    v = _diffs(df, model, kg); n = len(v)
    return v.mean(), st.t.ppf(0.975, n - 1) * v.std(ddof=1) / np.sqrt(n), (v > 0).mean()


# ════════════════════════════ 1 · KG lift by graph and model ════════════════════════════
def fig_kg_lift_bars(df):
    fig, ax = plt.subplots(figsize=(9.0, 5.4)); fig.subplots_adjust(left=0.09, right=0.97, top=0.84, bottom=0.12)
    ax.grid(axis='y', color='#ececec', lw=0.8); ax.set_axisbelow(True)
    w = 0.26; x = np.arange(len(KGS))
    for mi, m in enumerate(MODELS):
        off = (mi - 1) * w
        means = [_lift_ci(df, m, k)[0] for k in KGS]; cis = [_lift_ci(df, m, k)[1] for k in KGS]
        ax.bar(x + off, means, w, color=MCOL[m], label=m, edgecolor='white', lw=0.4,
               yerr=cis, error_kw=dict(elinewidth=1.2, capsize=3, ecolor='#333'))
    ax.set_xticks(x); ax.set_xticklabels([KLAB[k] for k in KGS], fontsize=12)
    ax.set_ylabel('Δ MRR vs no-KG  (95% CI)', fontsize=11); ax.set_ylim(0, None)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=3, fontsize=11, handlelength=1.2)
    ax.set_title('KG lift by graph and model', fontsize=14, fontweight='bold', pad=26)
    _footer(fig); return fig


# ════════════════════════════ 2 · per-disease lift box (faceted) ════════════════════════════
def fig_lift_box(df):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.8), sharex=True)
    fig.subplots_adjust(left=0.05, right=0.985, top=0.88, bottom=0.13, wspace=0.08)
    rng = np.random.default_rng(5)
    for ax, m in zip(axes, MODELS):
        ax.axvline(0, color='#c4c4c4', lw=0.9, ls='--', zorder=1)
        for j, kg in enumerate(KGS):
            y = 2 - j; v = _diffs(df, m, kg); mu = v.mean(); ph = (v > 0).mean() * 100
            ax.boxplot(v, positions=[y], vert=False, widths=0.5, patch_artist=True, showfliers=False, zorder=2,
                       medianprops=dict(color='#1a1a1a', lw=1.2),
                       boxprops=dict(facecolor=KCOL[kg], alpha=0.22, edgecolor=KCOL[kg], lw=1.1),
                       whiskerprops=dict(color=KCOL[kg], lw=1.1), capprops=dict(color=KCOL[kg], lw=1.1))
            ax.scatter(v, y + rng.uniform(-0.18, 0.18, len(v)), s=11, color=KCOL[kg], alpha=0.7,
                       edgecolor='white', linewidth=0.25, zorder=3)
            ax.text(0.69, y + 0.34, f'+{mu:.2f} · {ph:.0f}%', fontsize=9.5, color=KCOL[kg], fontweight='bold', ha='right')
        ax.set_yticks([2, 1, 0]); ax.set_yticklabels([KLAB[k] for k in KGS] if m == MODELS[0] else [], fontsize=10)
        ax.set_ylim(-0.55, 2.7); ax.set_xlim(-0.55, 0.72)
        ax.set_xlabel('Δ MRR vs no-KG  (dot = disease)', fontsize=9.5)
        ax.set_title(m, fontsize=12, fontweight='bold', pad=4)
    _footer(fig); return fig


# ════════════════════════════ 3 · rank distribution (faceted) ════════════════════════════
def _rbins(s):
    r = s['rank']
    return [(r == 1).mean(), (r == 2).mean(), (r == 3).mean(), ((r >= 4) & (r <= 8)).mean(), (r >= 9).mean()]

def fig_rank_distribution(df):
    arms = ['no_kg'] + KGS; alab = {'no_kg': 'No KG', **KLAB}
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 5.0), sharex=True)
    fig.subplots_adjust(left=0.05, right=0.985, top=0.80, bottom=0.11, wspace=0.08)
    for ax, m in zip(axes, MODELS):
        for i, a in enumerate(arms):
            vals = _rbins(df[(df.Model == m) & (df.arm == a)]); left = 0; y = len(arms) - 1 - i
            for v, c in zip(vals, RANK_COL):
                ax.barh(y, v, left=left, height=0.72, color=c, edgecolor='white', lw=0.6)
                if v > 0.05:
                    ax.text(left + v/2, y, f'{v*100:.0f}', ha='center', va='center', fontsize=8.5,
                            color='white' if c in ('#1F6F5E', '#C2705A') else '#333')
                left += v
        ax.set_yticks(range(len(arms))); ax.set_yticklabels([alab[a] for a in arms][::-1] if m == MODELS[0] else [], fontsize=10)
        ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_xlabel('Fraction of queries', fontsize=10)
        ax.set_title(m, fontsize=12, fontweight='bold', pad=4)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in RANK_COL]
    fig.legend(handles, RANK_BINS, ncol=5, loc='upper center', bbox_to_anchor=(0.5, 0.96), fontsize=10, handlelength=1.1)
    _footer(fig); return fig


# ════════════════════════════ 4 · which quality dimension tracks lift ════════════════════════════
def fig_dim_tracks_lift(df, summary_csv, min_range=0.02):
    b = pd.read_csv(summary_csv).set_index('KG')
    meanlift = {k: np.mean([_lift_ci(df, m, k)[0] for m in MODELS]) for k in KGS}
    order = sorted(KGS, key=lambda k: -meanlift[k]); y = np.array([meanlift[k] for k in KGS])
    rows = []
    for c in b.columns:
        x = b.loc[KGS, c].astype(float).values
        if np.isnan(x).any() or (x.max() - x.min()) < min_range: continue
        r = st.pearsonr(x, y)[0]
        match = sorted(KGS, key=lambda k: -b.loc[k, c]) == order
        rows.append((c, r, match))
    T = pd.DataFrame(rows, columns=['dim', 'r', 'match']).sort_values('r')
    fig, ax = plt.subplots(figsize=(7.2, 7.4)); fig.subplots_adjust(left=0.30, right=0.96, top=0.92, bottom=0.08)
    yy = np.arange(len(T))
    ax.barh(yy, T['r'], color=['#2E8B57' if m else '#c4c4c4' for m in T['match']], edgecolor='white', lw=0.5)
    ax.axvline(0, color='#666', lw=0.8); ax.set_yticks(yy); ax.set_yticklabels(T['dim'], fontsize=9)
    ax.set_xlim(-1.05, 1.05); ax.set_xlabel('Pearson r with mean lift  (n = 3 KGs)', fontsize=10)
    for t in ax.get_yticklabels():
        if t.get_text() == 'Coverage [dim]': t.set_fontweight('bold'); t.set_color('#1b5e34')
    h1 = plt.Rectangle((0, 0), 1, 1, color='#2E8B57'); h2 = plt.Rectangle((0, 0), 1, 1, color='#c4c4c4')
    ax.legend([h1, h2], ['KG order matches lift order', 'order differs'], loc='lower right', fontsize=8.5)
    ax.text(-1.0, len(T) - 0.6, 'exploratory: n = 3 KGs, no inference', fontsize=8, color='#999', style='italic')
    ax.set_title('Quality dimensions vs KG lift', fontsize=13.5, fontweight='bold', pad=10)
    return fig


# ════════════════════════════ 5 · lift vs coverage, faceted by LLM ════════════════════════════
def fig_lift_vs_coverage(df, summary_csv, xdim='Coverage [dim]'):
    b = pd.read_csv(summary_csv).set_index('KG'); cov = {k: float(b.loc[k, xdim]) for k in KGS}
    order = sorted(KGS, key=lambda k: cov[k])
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.80, bottom=0.15, wspace=0.1)
    for ax, m in zip(axes, MODELS):
        xs = np.array([cov[k] for k in order]); ys = np.array([_lift_ci(df, m, k)[0] for k in order])
        es = [_lift_ci(df, m, k)[1] for k in order]
        s, i = np.polyfit(xs, ys, 1); xl = np.array([xs.min(), xs.max()])
        ax.plot(xl, s*xl + i, color='#bdbdbd', lw=1.1, ls='--', zorder=2)
        ax.errorbar(xs, ys, yerr=es, fmt=VMK[m], ms=8, color=VCOL[m], mec='white', mew=0.7,
                    ecolor=VCOL[m], elinewidth=1.1, capsize=3, zorder=3)
        for k in KGS:
            ax.annotate(KLAB[k], (cov[k], _lift_ci(df, m, k)[0]), xytext=(5, -9), textcoords='offset points',
                        fontsize=8, color='#777')
        ax.set_xlim(0.44, 0.61); ax.set_xlabel('coverage', fontsize=10)
        ax.set_title(m, fontsize=12, fontweight='bold', color=VCOL[m], pad=6)
    axes[0].set_ylabel('Downstream lift (ΔMRR)', fontsize=10)
    fig.suptitle('Lift vs KG coverage, by LLM', fontsize=13, fontweight='bold', x=0.07, ha='left', y=0.97)
    return fig


# ════════════════════════════ 6 · evidence-quality (lollipop + calibration) ════════════════════════════
def fig_evidence_quality(df):
    kg = df[df.condition == 'kg']
    cats = [('consistent', 'Consistent\nevidence', '#2E8B57'),
            ('conflicting', 'Conflicting', '#E8A800'),
            ('insufficient', 'Insufficient /\nmissing', '#CB3B2E')]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.4), gridspec_kw=dict(width_ratios=[1.15, 1]))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.84, bottom=0.13, wspace=0.55)
    # left: lollipop of mean rank + side annotation
    for i, (key, lab, col) in enumerate(cats):
        s = kg[kg.pos_agreement == key]; y = len(cats) - 1 - i
        if len(s) == 0: continue
        mr = s['rank'].mean(); p1 = (s['rank'] == 1).mean() * 100; n = len(s)
        axL.plot([0, mr], [y, y], color=col, lw=4, solid_capstyle='round', zorder=2)
        axL.scatter(mr, y, s=520, color=col, zorder=3)
        axL.text(mr, y, f'{mr:.1f}', ha='center', va='center', color='white', fontsize=10, fontweight='bold')
        axL.annotate(f'{p1:.0f}% ranked #1\n(n={n:,})', xy=(8.6, y), ha='right', va='center',
                     fontsize=10.5, color=col, fontweight='bold')
    axL.set_yticks(range(len(cats))); axL.set_yticklabels([c[1] for c in cats][::-1], fontsize=10.5)
    axL.set_xlim(0, 11.5); axL.set_xticks(range(0, 9)); axL.set_ylim(-0.6, len(cats) - 0.4)
    axL.set_xlabel('Mean rank of the true drug   (1 = best, 8 = worst)', fontsize=10)
    axL.set_title('Mean rank of the true drug by KG evidence', fontsize=11.5,
                  fontweight='bold', loc='left', color='#2b3a55', pad=6)
    # right: confidence calibration
    g = kg.groupby('pos_confidence')['hit@1'].mean()
    x = g.index.values; yv = g.values
    axR.grid(axis='y', color='#ececec', lw=0.8); axR.set_axisbelow(True)
    axR.plot(x, yv, color='#9a9a9a', lw=1.6, zorder=2)
    axR.scatter(x, yv, s=300, c=[plt.cm.RdYlGn(v) for v in yv], edgecolor='white', linewidth=1.3, zorder=3)
    for xi, yi in zip(x, yv):
        axR.annotate(f'{yi*100:.0f}%', (xi, yi), xytext=(0, 12), textcoords='offset points', ha='center',
                     fontsize=9.5, fontweight='bold', color='#333')
    axR.set_xticks([1, 2, 3, 4, 5]); axR.set_xlim(0.6, 5.4); axR.set_ylim(-0.04, 1.1)
    axR.set_xlabel("Model's self-reported confidence (1–5)", fontsize=10); axR.set_ylabel('Hits@1', fontsize=10)
    axR.set_title('Hits@1 by model confidence', fontsize=11.5, fontweight='bold', loc='left', pad=6)
    fig.suptitle('True-drug rank and Hits@1 vs KG evidence', fontsize=14, fontweight='bold', y=0.97)
    _footer(fig, '116-set · 3 seeds × 2 shuffles · GPT + Gemini + Llama, KG arms · model self-reports')
    return fig


def regenerate(runs_dir, summary_csv, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    df = load_runs(runs_dir)
    figs = {'fig1_kg_lift_bars': fig_kg_lift_bars(df),
            'fig2_lift_box': fig_lift_box(df),
            'fig3_rank_distribution': fig_rank_distribution(df),
            'fig4_dim_tracks_lift': fig_dim_tracks_lift(df, summary_csv),
            'fig5_lift_vs_coverage': fig_lift_vs_coverage(df, summary_csv),
            'fig6_evidence_quality': fig_evidence_quality(df)}
    for name, fig in figs.items():
        fig.savefig(os.path.join(out_dir, name + '.png'), bbox_inches='tight', facecolor='white')
        fig.savefig(os.path.join(out_dir, name + '.pdf'), bbox_inches='tight', facecolor='white')
    return figs
