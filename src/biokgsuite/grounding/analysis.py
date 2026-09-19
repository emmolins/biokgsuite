"""Statistics for the grounding experiment: one tidy run frame in, tidy result
tables out.

Every function takes the run frame written by ``runner`` (one row per query x
arm x model x seed x shuffle, carrying ``pair_id``, ``family`` and ``arm``) and
returns a DataFrame that a notebook writes through ``Outputs`` and a figure in
``biokgsuite.viz.grounding`` draws from. Nothing here reads files except
``load_runs``, and nothing here knows a model name, a graph name or a pair
count: those arrive as arguments.

The experimental unit is the drug-disease pair. Reciprocal ranks are averaged
to one value per (model, pair, arm) across seeds and shuffles before any
interval is computed, and the bootstrap resamples pairs, because rows within a
pair are not independent.
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# The key that makes a run row unique. A duplicate means a rerun was appended
# rather than replaced, and every mean downstream would be silently weighted.
RUN_KEY = ['family', 'model', 'arm', 'seed', 'shuffle', 'pair_id']

# Display name for the baseline arm when the caller's label map omits it. The
# arm name itself ('none') is part of the run-row vocabulary, not a display
# choice, so this is the one label the module supplies.
NONE_LABEL = 'No KG'

# Prefix for the per-model MRR columns in the difficulty table. Model tags
# contain ':' and '.', so they are namespaced rather than used bare, and the
# prefix is what lets the figures find them without a list of models.
MODEL_COL_PREFIX = 'mrr:'


def _label(labels, key):
    if key == 'none':
        return labels.get('none', NONE_LABEL)
    return labels.get(key, key)


def load_runs(run_dir, *, families=None, models=None, allow_mock=False):
    """Concatenate the per-(family, model, seed) run CSVs under ``run_dir``.

    Filtering is on the ``family`` and ``model`` columns, not on file names,
    so a file whose slug and content disagree cannot slip through. Rows the
    runner produced in mock mode (a deterministic stand-in for a model, used
    to exercise the pipeline) are refused unless ``allow_mock`` is set, so a
    rehearsal can never be analysed as a result.
    """
    files = sorted(Path(run_dir).glob('*.csv'))
    if not files:
        raise FileNotFoundError(f'no run files under {run_dir}')
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    if 'mode' in df.columns and not allow_mock and (df['mode'] == 'mock').any():
        bad = sorted(df.loc[df['mode'] == 'mock', 'model'].unique())
        raise ValueError(f'{run_dir} holds mock-mode rows for {bad}; these are rehearsal '
                         f'output, not results. Delete them or pass allow_mock=True')
    if families is not None:
        df = df[df['family'].isin(list(families))]
    if models is not None:
        df = df[df['model'].isin(list(models))]
    return df.reset_index(drop=True)


def validate_runs(df):
    """Assert that the frame can be analysed as a paired design.

    Two checks. No (family, model, arm, seed, shuffle, pair_id) may repeat, and
    within each (family, model, seed, shuffle) every arm must cover exactly the
    same set of pairs, otherwise the per-pair differences between arms would be
    computed on different pairs. Returns the frame unchanged so it can sit in a
    call chain.
    """
    missing = [c for c in RUN_KEY + ['reciprocal_rank'] if c not in df.columns]
    if missing:
        raise ValueError(f'run frame lacks columns {missing}')
    dup = df.duplicated(RUN_KEY)
    if dup.any():
        first = df.loc[dup, RUN_KEY].iloc[0].to_dict()
        raise ValueError(f'{int(dup.sum())} duplicated run rows, first {first}')
    for block, grp in df.groupby(['family', 'model', 'seed', 'shuffle'], sort=False):
        pairs_by_arm = grp.groupby('arm')['pair_id'].apply(lambda s: frozenset(s))
        if pairs_by_arm.nunique() > 1:
            union = frozenset().union(*pairs_by_arm)
            short = {arm: sorted(union - p) for arm, p in pairs_by_arm.items() if p != union}
            raise ValueError(f'arms do not cover the same pairs in block {block}: '
                             f'missing {short}')
    return df


def chance_mrr(pool_size):
    """Expected reciprocal rank of a uniformly random ranking of ``pool_size``."""
    return float(np.mean([1.0 / k for k in range(1, pool_size + 1)]))


def diagnostics(df, pool_size):
    """Parse rate, failed calls and a miss-corrected MRR per (family, model, arm).

    An unparsed response is scored as a miss (rank = pool_size + 1), which
    deflates MRR. A miss costs more in an arm that would otherwise score high,
    so failures bias against the best arms and the reported lifts are
    conservative. ``mrr_corrected`` is what MRR would be if the failed rows had
    scored at the arm's own mean, which bounds the size of that bias.
    """
    # A recorded error is a call that never returned a ranking; a parse failure
    # is a call that returned one the parser could not read. Both are misses.
    if 'error' in df.columns:
        failed_call = df['error'].notna() & (df['error'].astype(str).str.strip() != '')
    else:
        failed_call = pd.Series(False, index=df.index)
    diag = (df.assign(_failed_call=failed_call)
              .groupby(['family', 'model', 'arm'], sort=False)
              .agg(n=('reciprocal_rank', 'size'),
                   parse_rate=('parsed', 'mean'),
                   mrr=('reciprocal_rank', 'mean'),
                   n_error=('_failed_call', 'sum'))
              .reset_index())
    diag['n_error'] = diag['n_error'].astype(int)
    diag['n_failed'] = ((1 - diag['parse_rate']) * diag['n']).round().astype(int)
    miss_rr = 1.0 / (pool_size + 1)
    diag['mrr_corrected'] = diag['mrr'] + (diag['n_failed'] / diag['n']) * (diag['mrr'] - miss_rr)
    diag['shift'] = diag['mrr_corrected'] - diag['mrr']
    return diag


def per_query(df):
    """One row per (model, pair_id), one column per arm, mean reciprocal rank
    over seeds and shuffles. Pairs missing from any arm are dropped, so every
    arm is scored on the same pairs."""
    return (df.pivot_table(index=['model', 'pair_id'], columns='arm',
                           values='reciprocal_rank', aggfunc='mean')
              .dropna())


def boot_ci(values, n_boot, rng):
    """Percentile 95% interval for the mean of ``values`` by resampling with
    replacement. Fewer than two values give (nan, nan) rather than a degenerate
    interval, so a thin subset shows up as missing instead of certain."""
    v = np.asarray(values, dtype=float)
    if len(v) < 2:
        return (np.nan, np.nan)
    idx = rng.integers(0, len(v), size=(n_boot, len(v)))
    means = v[idx].mean(axis=1)
    return (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


def _arm_rows(pq, arms, n_boot, rng, labels):
    rows = []
    for model in pq.index.get_level_values('model').unique():
        p = pq.xs(model, level='model')
        for arm in [a for a in arms if a in p.columns]:
            lo, hi = boot_ci(p[arm], n_boot, rng)
            row = dict(model=model, model_label=_label(labels, model),
                       arm=arm, arm_label=_label(labels, arm),
                       n_units=len(p), mrr=float(p[arm].mean()), lo=lo, hi=hi,
                       lift=np.nan, lift_lo=np.nan, lift_hi=np.nan)
            if arm != 'none' and 'none' in p.columns:
                d = (p[arm] - p['none']).to_numpy()
                dlo, dhi = boot_ci(d, n_boot, rng)
                row.update(lift=float(d.mean()), lift_lo=dlo, lift_hi=dhi)
            rows.append(row)
    return rows


def arm_table(df, family, arms, n_boot, rng, labels):
    """MRR with bootstrap interval per (model, arm) within one family, plus the
    paired lift over the 'none' arm for every KG arm."""
    sub = df[df['family'] == family]
    if sub.empty:
        raise ValueError(f'no rows for family {family!r}')
    out = pd.DataFrame(_arm_rows(per_query(sub), arms, n_boot, rng, labels))
    out.insert(0, 'family', family)
    return out


def contrast_table(df, family, kg_pairs, n_boot, rng, labels):
    """Paired difference in per-pair MRR between two KG arms, per model.

    ``separates`` is True when the interval excludes zero. The difference is
    taken on the same pairs for both arms, which is what makes it a contrast
    between graphs rather than between subsets.
    """
    sub = df[df['family'] == family]
    if sub.empty:
        raise ValueError(f'no rows for family {family!r}')
    pq, rows = per_query(sub), []
    for model in pq.index.get_level_values('model').unique():
        p = pq.xs(model, level='model')
        for a, b in kg_pairs:
            if a not in p.columns or b not in p.columns:
                continue
            d = (p[a] - p[b]).to_numpy()
            lo, hi = boot_ci(d, n_boot, rng)
            rows.append(dict(family=family, model=model, model_label=_label(labels, model),
                             contrast=f'{_label(labels, a)} - {_label(labels, b)}',
                             kg_a=a, kg_b=b, n_units=len(p), diff=float(d.mean()),
                             lo=lo, hi=hi, separates=bool(lo > 0 or hi < 0)))
    return pd.DataFrame(rows)


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q), same order as the input."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return np.empty(0)
    order = np.argsort(p)
    q = np.empty(n)
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, p[order[i]] * n / (i + 1))
        q[order[i]] = running
    return q


def _wilcoxon_p(d):
    # Wilcoxon discards zero differences; when every difference is zero the
    # test has no data, and the honest p-value for "no difference" is 1.
    if not np.any(d != 0):
        return 1.0
    return float(stats.wilcoxon(d, alternative='two-sided').pvalue)


def significance_tests(df_by_subset, kg_pairs, labels):
    """Wilcoxon signed-rank and paired t-test for every KG contrast, per model
    and subset, with Benjamini-Hochberg FDR applied within each subset.

    Wilcoxon is the headline test because per-pair differences are bounded and
    skewed; the t-test is reported beside it so the two can be compared. FDR is
    applied within a subset because that is the family of tests a reader would
    scan together.
    """
    rows = []
    for subset, frame in df_by_subset.items():
        pq = per_query(frame)
        for model in pq.index.get_level_values('model').unique():
            p = pq.xs(model, level='model')
            for a, b in kg_pairs:
                if a not in p.columns or b not in p.columns:
                    continue
                d = (p[a] - p[b]).to_numpy()
                rows.append(dict(
                    subset=subset, model=model, model_label=_label(labels, model),
                    contrast=f'{_label(labels, a)} - {_label(labels, b)}',
                    kg_a=a, kg_b=b, n=len(d), diff=float(d.mean()),
                    p_wilcoxon=_wilcoxon_p(d),
                    p_ttest=float(stats.ttest_rel(p[a], p[b]).pvalue)))
    tests = pd.DataFrame(rows)
    if tests.empty:
        return tests
    tests['q_bh'] = tests.groupby('subset')['p_wilcoxon'].transform(bh_fdr)
    tests['sig_raw'] = tests['p_wilcoxon'] < 0.05
    tests['sig_fdr'] = tests['q_bh'] < 0.05
    return tests


def capacity_table(df, params_b, arms, n_boot, rng, labels):
    """Arm table for a model ladder, sorted by parameter count.

    ``params_b`` maps model tag to size in billions and also selects which
    models belong to the ladder, so a frame holding other models is fine.
    """
    sub = df[df['model'].isin(list(params_b))]
    if sub.empty:
        raise ValueError('no rows for any model in params_b')
    cap = pd.DataFrame(_arm_rows(per_query(sub), arms, n_boot, rng, labels))
    cap['params_b'] = cap['model'].map(params_b).astype(float)
    arm_rank = {a: i for i, a in enumerate(arms)}
    cap['_arm_rank'] = cap['arm'].map(arm_rank)
    cap = cap.sort_values(['params_b', '_arm_rank']).drop(columns='_arm_rank')
    return cap.reset_index(drop=True)


def scaling_table(df, params_b, grounded_pairs_by_kg, n_boot, rng):
    """Per rung of a model ladder: no-KG MRR, KG MRR pooled over every KG arm,
    KG MRR on grounded arms only, and the paired lift of each over no-KG.

    Splitting pooled from grounded is what separates two readings of a flat
    lift: a graph with no evidence for the pair cannot help however capable
    the model, so the grounded-arm lift is the one that speaks to capability.
    ``grounded_pairs_by_kg`` maps an arm to the pair ids whose positive has a
    dossier in that graph (``grounded_<kg>`` in the benchmark frame).
    """
    kg_rows = df[(df['arm'] != 'none') & df['model'].isin(list(params_b))].copy()
    kg_rows['grounded'] = [pid in grounded_pairs_by_kg.get(arm, set())
                           for pid, arm in zip(kg_rows['pair_id'], kg_rows['arm'])]
    rows = []
    for model in sorted(params_b, key=params_b.get):
        sub = df[df['model'] == model]
        if sub.empty:
            continue
        none = sub[sub['arm'] == 'none'].groupby('pair_id')['reciprocal_rank'].mean()
        kg = kg_rows[kg_rows['model'] == model]
        pooled = kg.groupby('pair_id')['reciprocal_rank'].mean()
        grd = kg[kg['grounded']].groupby('pair_id')['reciprocal_rank'].mean()
        lift_pooled = (pooled - none).dropna().to_numpy()
        lift_grd = (grd - none).dropna().to_numpy()
        lo_p, hi_p = boot_ci(lift_pooled, n_boot, rng)
        lo_g, hi_g = boot_ci(lift_grd, n_boot, rng)
        rows.append(dict(
            model=model, params_b=float(params_b[model]),
            mrr_none=float(none.mean()), mrr_kg=float(pooled.mean()),
            mrr_kg_grounded=float(grd.mean()) if len(grd) else np.nan,
            lift_pooled=float(lift_pooled.mean()) if len(lift_pooled) else np.nan,
            lift_pooled_lo=lo_p, lift_pooled_hi=hi_p,
            lift_grounded=float(lift_grd.mean()) if len(lift_grd) else np.nan,
            lift_grounded_lo=lo_g, lift_grounded_hi=hi_g,
            n_pairs=len(lift_pooled), n_pairs_grounded=len(lift_grd)))
    return pd.DataFrame(rows).sort_values('params_b').reset_index(drop=True)


def absence_penalty_table(df, n_boot, rng, labels):
    """What a graph does to a positive it has nothing to say about.

    In a KG arm every candidate the graph holds carries a dossier; a positive
    the graph holds no gene or pathway for reaches the model with none while
    its distractors carry theirs. Per (model, arm) this splits the pairs by
    whether the positive's dossier was shown (``pos_has_evidence``, recorded
    on the row, so this reads what the model saw rather than a benchmark
    flag) and reports the paired lift over the ``none`` arm for each half,
    with bootstrap intervals over pairs. The absent-positive lift is the
    penalty a coverage gap imposes at deployment; a negative value means the
    graph ranked its own blind spots below chance. Pairs whose eight
    candidates were all blank are left out of the absent half, since there
    the arm reduces to ``none``.
    """
    none = (df[df['arm'] == 'none'].groupby(['model', 'pair_id'])['reciprocal_rank'].mean()
              .rename('none'))
    kg = df[df['arm'] != 'none'].copy()
    kg['pos_has_evidence'] = kg['pos_has_evidence'].astype(bool)
    kg['any_evidence'] = kg['n_with_evidence'].astype(int) > 0
    rows = []
    for (model, arm), sub in kg.groupby(['model', 'arm'], sort=False):
        pq = sub.groupby('pair_id').agg(rr=('reciprocal_rank', 'mean'),
                                        shown=('pos_has_evidence', 'all'),
                                        any_evidence=('any_evidence', 'any'))
        pq = pq.join(none.xs(model, level='model'), how='inner')
        halves = {'shown': pq[pq['shown']],
                  'absent': pq[~pq['shown'] & pq['any_evidence']]}
        for half, h in halves.items():
            d = (h['rr'] - h['none']).to_numpy()
            lo, hi = boot_ci(d, n_boot, rng)
            rows.append(dict(model=model, model_label=_label(labels, model),
                             arm=arm, arm_label=_label(labels, arm), positive=half,
                             n_pairs=len(h), mrr=float(h['rr'].mean()) if len(h) else np.nan,
                             mrr_none=float(h['none'].mean()) if len(h) else np.nan,
                             lift=float(d.mean()) if len(d) else np.nan, lift_lo=lo, lift_hi=hi))
    return pd.DataFrame(rows)


BENCHMARK_COVARIATES = ['drugbank_id', 'mondo_id', 'drug_name', 'disease_name',
                        'area', 'modality', 'approval_year']


def difficulty_table(df, benchmark):
    """One row per pair: how hard it is to rank, pooled over models, seeds and
    shuffles on the KG arms, joined to the benchmark covariates.

    ``mrr`` pools every KG arm, including arms whose graph has no evidence
    for the pair; ``mrr_grounded`` uses only arms whose graph does. The gap
    between the two is the empty-dossier effect, and both are needed to tell
    "hard because the graphs have nothing on it" from "hard even with
    evidence". ``n_kg_covered`` and ``n_kg_grounded`` come from the benchmark.
    """
    kg = df[df['arm'] != 'none'].copy()
    kg_arms = sorted(kg['arm'].unique())
    grd_cols = {arm: f'grounded_{arm}' for arm in kg_arms}
    missing = [c for c in grd_cols.values() if c not in benchmark.columns]
    if missing:
        raise KeyError(f'benchmark frame lacks {missing}; rebuild the pair status and rerun 01_benchmark')
    bench = benchmark.drop_duplicates('pair_id').set_index('pair_id')

    grounded = pd.Series(False, index=kg.index)
    for arm, col in grd_cols.items():
        hit = kg['arm'] == arm
        grounded[hit] = kg.loc[hit, 'pair_id'].map(bench[col]).fillna(False).astype(bool).to_numpy()
    kg['grounded'] = grounded.to_numpy()

    per = (kg.groupby('pair_id')
             .agg(mrr=('reciprocal_rank', 'mean'), hit1=('hit1', 'mean'),
                  pos_confidence=('pos_confidence', 'mean'),
                  n_runs=('reciprocal_rank', 'size'))
             .reset_index())
    per['frac_consistent'] = (per['pair_id']
                              .map((kg['pos_agreement'] == 'consistent').groupby(kg['pair_id']).mean()))
    per['mrr_grounded'] = per['pair_id'].map(
        kg[kg['grounded']].groupby('pair_id')['reciprocal_rank'].mean())
    per['mrr_none'] = per['pair_id'].map(
        df[df['arm'] == 'none'].groupby('pair_id')['reciprocal_rank'].mean())
    per['n_kg_grounded'] = per['pair_id'].map(bench[list(grd_cols.values())].sum(axis=1))
    per['n_kg_covered'] = (per['pair_id'].map(bench['n_kg_covered'])
                           if 'n_kg_covered' in bench.columns else np.nan)

    model_cols = []
    for model in kg['model'].unique():
        col = f'{MODEL_COL_PREFIX}{model}'
        per[col] = per['pair_id'].map(
            kg[kg['model'] == model].groupby('pair_id')['reciprocal_rank'].mean())
        model_cols.append(col)
    per['model_spread'] = per[model_cols].max(axis=1) - per[model_cols].min(axis=1)

    keep = [c for c in BENCHMARK_COVARIATES if c in bench.columns]
    per = per.merge(bench[keep], left_on='pair_id', right_index=True, how='left')
    return per.sort_values('mrr').reset_index(drop=True)


def model_columns(per):
    """The per-model MRR columns of a difficulty table, as {column: model}."""
    return {c: c[len(MODEL_COL_PREFIX):] for c in per.columns if c.startswith(MODEL_COL_PREFIX)}


def cross_model_corr(per):
    """Pearson correlation of per-pair MRR between models. Off-diagonal values
    near one mean difficulty is a property of the pair, not the model."""
    cols = model_columns(per)
    corr = per[list(cols)].corr()
    return corr.rename(index=cols, columns=cols)


def mean_pairwise_corr(corr):
    """Mean of the upper triangle of a correlation matrix, the single number
    the intrinsic-difficulty figure quotes."""
    vals = corr.to_numpy()
    if len(vals) < 2:
        return np.nan
    return float(np.mean([vals[i, j] for i, j in combinations(range(len(vals)), 2)]))
