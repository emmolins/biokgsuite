"""Selection-bias audit: does the kept subset of a benchmark differ from the
rows that were dropped?

Dataset-agnostic. The caller supplies a frame with a boolean kept column, the
categorical covariates to test level by level, and the continuous covariates
to compare by rank test. The figures that draw these results live in
``biokgsuite.viz.grounding``.
"""
import numpy as np
import pandas as pd
from scipy import stats


def odds_ratio(mask, kept):
    """Odds of being kept given ``mask`` is True, as (or, lo, hi, p).

    The interval uses the Haldane correction (0.5 added to every cell) so a
    level that was never dropped, or never kept, still yields a finite
    estimate rather than an infinite one. The p-value is Fisher's exact test on
    the uncorrected counts, because the correction is for estimation, not for
    inference.
    """
    mask = np.asarray(mask, dtype=bool)
    kept = np.asarray(kept, dtype=bool)
    a = int((mask & kept).sum())
    b = int((mask & ~kept).sum())
    c = int((~mask & kept).sum())
    d = int((~mask & ~kept).sum())
    ah, bh, ch, dh = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    lor = np.log(ah * dh / (bh * ch))
    se = np.sqrt(1 / ah + 1 / bh + 1 / ch + 1 / dh)
    p = float(stats.fisher_exact([[a, b], [c, d]])[1])
    return (float(np.exp(lor)), float(np.exp(lor - 1.96 * se)),
            float(np.exp(lor + 1.96 * se)), p)


def run_bias_tests(df, kept_col, cat_vars, cont_vars, min_n=8):
    """Tidy results table, one row per categorical level (one-vs-rest odds
    ratio) and one per continuous covariate (Mann-Whitney on kept vs dropped).

    ``effect`` is the odds ratio for a categorical row and the difference of
    medians (kept minus dropped) for a continuous row; the two medians are
    also carried so the sign can be read back. Levels with fewer than
    ``min_n`` rows are skipped because their odds ratios are dominated by the
    correction rather than by the data.
    """
    kept = df[kept_col].astype(bool)
    rows = []
    for var in cat_vars:
        for level, n in df[var].value_counts().items():
            if n < min_n:
                continue
            orr, lo, hi, p = odds_ratio(df[var] == level, kept)
            rows.append(dict(covariate=var, trait=f'{var} = {level}', kind='categorical',
                             n=int(n), effect=orr, ci_low=lo, ci_high=hi, p=p,
                             median_kept=np.nan, median_dropped=np.nan))
    for var in cont_vars:
        a = df.loc[kept, var].dropna()
        b = df.loc[~kept, var].dropna()
        if len(a) == 0 or len(b) == 0:
            continue
        p = float(stats.mannwhitneyu(a, b).pvalue)
        rows.append(dict(covariate=var, trait=f'{var} (kept vs dropped)', kind='continuous',
                         n=int(len(a) + len(b)), effect=float(a.median() - b.median()),
                         ci_low=np.nan, ci_high=np.nan, p=p,
                         median_kept=float(a.median()), median_dropped=float(b.median())))
    return pd.DataFrame(rows, columns=['covariate', 'trait', 'kind', 'n', 'effect',
                                       'ci_low', 'ci_high', 'p', 'median_kept',
                                       'median_dropped'])
