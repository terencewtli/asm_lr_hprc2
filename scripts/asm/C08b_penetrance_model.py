#!/usr/bin/env python3
"""C08b_penetrance_model.py -- replace the donor-tier split with a per-donor propensity offset,
so penetrance is estimated on ALL donors, continuously, with a confidence interval.

Why (RESULTS §7, §10; PROGRESS "ASM replication outputs"):
  C07d estimates penetrance as n_asm_rep / n_tested_rep over the 105 replication-tier donors and
  calls a region "replicating" at binomial p < 1e-3 against a single scalar P0 = 0.0454. Three
  problems, all fixed here:
    1. the scalar P0 gives a hard detection floor of 13/105 donors = 12.4% penetrance;
    2. 96 donors (69 discovery + 22 inflated + 5 excluded) contribute nothing to penetrance, and
       because lambda and yield correlate at rho = 0.974 the 69 discarded discovery donors are
       exactly the cleanest ones;
    3. P0 is estimated over all candidates, so true signal inflates the null.

Model. For donor i and region j, with b_i the donor's propensity to call ASM at a null region:
    logit P(call_ij) = theta_j + logit(b_i)
theta_j is the region's excess on the log-odds scale, adjusted for how call-happy each donor is.
b_i absorbs lambda, domain depth, chemistry and coverage without having to model them explicitly
-- it is measured, not predicted (the donor-covariate regression is reported as a diagnostic).

b_i is estimated on a NULL SET of regions rather than assumed. The default null set is ALL
candidates, which is deliberately conservative (it contains the true signal, so b_i is biased up
and the test is under-powered) and is the per-donor analogue of C07d's scalar P0 = 0.0454.
`--null-quantile q` instead estimates b_i on regions whose raw penetrance is at or below the q-th
quantile.

**Do not iterate this to convergence.** An earlier version re-estimated b_i on the regions left
after dropping q < 0.05 and repeated: the null set eroded 119,024 -> 73,832 over four passes with
b_i still falling (0.0442 -> 0.0145) and 37.9% of candidates called. Each pass removes regions,
which lowers b_i, which makes more regions significant. There is no fixed point. The number of
called regions is strongly sensitive to this choice (13.9% / 29.6% / 43.6% at null quantile
1.00 / 0.75 / 0.50), so the sensitivity table is part of the output, not a footnote.

Note that genomic control in C07b zeroes out the most inflated donors entirely -- HG02040,
HG01150, HG02129, HG03540 and NA18608 call nothing at any candidate, so b_i = 0 and they
contribute no information. The effective donor count is ~196, not 201.

Reported per region:
  penetrance_raw    n_asm / n_tested over all donors (no adjustment)
  penetrance_adj    fitted P(call) at the median-propensity donor -- "penetrance in a typical
                    donor", the quantity to plot
  penetrance_lo/hi  95% CI on penetrance_adj (Wald on theta)
  theta, se_theta, p_excess (one-sided), q_excess (BH)

Outputs (results/asm/model/):
  penetrance_adjusted.tsv.gz     one row per candidate region
  donor_propensity.tsv           b_i per donor + covariates + the diagnostic regression
  penetrance_model_summary.tsv   comparison against the C07d tier estimate

Usage: C08b_penetrance_model.py [--min-tested 20] [--iters 4]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.stats import norm

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
REP = PROJDIR / 'results' / 'asm' / 'replication'
QC = PROJDIR / 'results' / 'qc' / 'data'
OUT = PROJDIR / 'results' / 'asm' / 'model'
EPS = 1e-9


def bh(p: np.ndarray) -> np.ndarray:
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pp = p[ok]
    o = np.argsort(pp)
    n = len(pp)
    qq = np.empty(n)
    qq[o] = np.minimum.accumulate((pp[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    q[ok] = np.minimum(qq, 1.0)
    return q


def fit_theta(K: np.ndarray, M: np.ndarray, off: np.ndarray, lo=-12.0, hi=12.0, it=50):
    """Vectorised bisection for theta_j solving  K_j = sum_{i tested} expit(theta_j + off_i)."""
    n = M.shape[0]
    a = np.full(n, lo)
    b = np.full(n, hi)
    for _ in range(it):
        mid = 0.5 * (a + b)
        s = (expit(mid[:, None] + off[None, :]) * M).sum(1)
        too_low = s < K
        a = np.where(too_low, mid, a)
        b = np.where(too_low, b, mid)
    th = 0.5 * (a + b)
    p = expit(th[:, None] + off[None, :]) * M
    info = (p * (1 - p) * M).sum(1)
    se = 1.0 / np.sqrt(np.maximum(info, EPS))
    return th, se


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-tested', type=int, default=20)
    ap.add_argument('--null-quantile', type=float, default=1.0,
                    help='estimate b_i on regions with raw penetrance <= this quantile (1.0 = all '
                         'candidates, the conservative default)')
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    tiers = pd.read_csv(REP / 'donor_tiers.tsv', sep='\t')
    files = {s: REP / 'calls' / ('%s.tsv.gz' % s) for s in tiers['sample']}
    samples = [s for s, f in files.items() if f.exists()]
    print('donors with re-test calls: %d (C07d used the 105 replication-tier only)' % len(samples), flush=True)

    calls = pd.concat([pd.read_csv(files[s], sep='\t', usecols=['sample', 'region_id', 'asm_rep'])
                       for s in samples], ignore_index=True)
    regions = pd.Index(sorted(calls.region_id.unique()))
    donors = pd.Index(samples)
    ri = regions.get_indexer(calls.region_id)
    di = donors.get_indexer(calls['sample'])
    M = np.zeros((len(regions), len(donors)), dtype=np.float32)
    A = np.zeros_like(M)
    M[ri, di] = 1.0
    A[ri, di] = calls.asm_rep.values.astype(np.float32)
    del calls
    K = A.sum(1)
    N = M.sum(1)
    print('matrix: %d regions x %d donors; median tested %d' % (len(regions), len(donors), np.median(N)), flush=True)

    raw = np.where(N > 0, K / np.maximum(N, 1), np.nan)

    def propensity(q):
        msk = np.ones(len(regions), bool) if q >= 1.0 else (raw <= np.nanquantile(raw, q))
        return (A[msk].sum(0) + 0.5) / (M[msk].sum(0) + 1.0), int(msk.sum())

    def fit_and_call(b):
        off = logit(np.clip(b, 1e-6, 1 - 1e-6))
        th, se = fit_theta(K, M, off)
        p = norm.sf(th / se)
        return th, se, p, bh(p)

    # sensitivity across null-set choices, reported alongside the headline
    sens = []
    for q in (1.0, 0.75, 0.5):
        b_q, nq = propensity(q)
        th_q, se_q, p_q, q_q = fit_and_call(b_q)
        sig = q_q < 0.05
        pen_q = expit(th_q + logit(np.median(b_q)))
        sens.append({'null_quantile': q, 'n_null_regions': nq, 'b_i_median': float(np.median(b_q)),
                     'n_called': int(sig.sum()), 'frac_called': float(sig.mean()),
                     'min_penetrance_called': float(np.nanmin(np.where(sig, pen_q, np.nan)))})
        print('  null q=%.2f: %d regions, b_i median %.4f -> %d called (%.2f%%), floor %.4f'
              % (q, nq, np.median(b_q), sig.sum(), 100 * sig.mean(), sens[-1]['min_penetrance_called']),
              flush=True)
    pd.DataFrame(sens).to_csv(OUT / 'penetrance_null_sensitivity.tsv', sep='\t', index=False,
                              float_format='%.5g')

    b_i, n_null = propensity(a.null_quantile)
    th, se, p, q = fit_and_call(b_i)
    zero = donors[b_i <= 0.002]
    print('\nheadline null-quantile=%.2f (%d regions); donors with b_i ~ 0 (GC zeroed them): %s'
          % (a.null_quantile, n_null, list(zero)), flush=True)

    b_ref = float(np.median(b_i))
    o_ref = logit(b_ref)
    pen_adj = expit(th + o_ref)
    pen_lo = expit(th - 1.96 * se + o_ref)
    pen_hi = expit(th + 1.96 * se + o_ref)

    out = pd.DataFrame({
        'region_id': regions, 'n_tested_all': N.astype(int), 'n_asm_all': K.astype(int),
        'penetrance_raw': np.where(N > 0, K / np.maximum(N, 1), np.nan),
        'theta': th, 'se_theta': se, 'penetrance_adj': pen_adj,
        'penetrance_lo': pen_lo, 'penetrance_hi': pen_hi,
        'p_excess': p, 'q_excess': q})
    out.loc[out.n_tested_all < a.min_tested, ['penetrance_adj', 'penetrance_lo', 'penetrance_hi',
                                              'p_excess', 'q_excess']] = np.nan
    out.to_csv(OUT / 'penetrance_adjusted.tsv.gz', sep='\t', index=False, float_format='%.5g')

    # donor propensity table + diagnostic: what does b_i track?
    dp = pd.DataFrame({'sample': donors, 'b_i': b_i, 'informative': b_i > 0.002,
                       'n_tested': M.sum(0).astype(int),
                       'n_asm': A.sum(0).astype(int)}).merge(tiers, on='sample', how='left')
    qc = pd.read_csv(QC / 'qc17' / 'donor_metrics_covariates.tsv', sep='\t')
    dp = dp.merge(qc[['sample', 'depth_constitutive', 'sex', 'coverage_ont']], on='sample', how='left')
    nl = pd.read_csv(QC / 'asm_null_vs_real_chr20.tsv', sep='\t')[['sample', 'read_sd', 'lambda_null']]
    dp = dp.merge(nl, on='sample', how='left')
    dp.to_csv(OUT / 'donor_propensity.tsv', sep='\t', index=False, float_format='%.5g')

    import statsmodels.formula.api as smf
    dd = dp[dp.informative].dropna(subset=['depth_constitutive', 'read_sd', 'lambda_gc']).copy()
    dd['logb'] = np.log(dd.b_i)
    diag = []
    for f in ['logb ~ np.log(lambda_gc)', 'logb ~ read_sd + depth_constitutive',
              'logb ~ np.log(lambda_gc) + read_sd + depth_constitutive + C(chemistry)']:
        r = smf.ols(f, data=dd).fit()
        diag.append((f, r.rsquared))
    print('\ndonor propensity b_i: median %.4f, range %.4f-%.4f (%.1fx spread)' % (
        b_ref, b_i.min(), b_i.max(), b_i.max() / b_i.min()))
    for f, r2 in diag:
        print('  R2=%.3f  %s' % (r2, f))

    # comparison against C07d
    c7 = pd.read_csv(REP / 'candidates_replication.tsv.gz', sep='\t',
                     usecols=['region_id', 'set', 'rep_class', 'penetrance_rep', 'n_tested_rep',
                              'zink_10kb', 'p_replicate'], low_memory=False)
    cmp = c7.merge(out, on='region_id', how='left')
    sig_new = cmp.q_excess < 0.05
    sig_old = cmp.p_replicate < 1e-3
    print('\n=== C07d tier estimate vs C08b adjusted ===')
    print('  regions called (old, p_replicate<1e-3):     %7d (%.2f%%)' % (sig_old.sum(), 100 * sig_old.mean()))
    print('  regions called (new, q_excess<0.05):        %7d (%.2f%%)' % (sig_new.sum(), 100 * sig_new.mean()))
    print('  gained %d, lost %d, shared %d' % ((sig_new & ~sig_old).sum(), (~sig_new & sig_old).sum(),
                                               (sig_new & sig_old).sum()))
    print('  median n_tested: old %.0f -> new %.0f  (+%.0f donors per region)' % (
        cmp.n_tested_rep.median(), cmp.n_tested_all.median(), cmp.n_tested_all.median() - cmp.n_tested_rep.median()))
    zk = cmp[cmp.set == 'zink']
    print('  Zink DMRs recovered: old %d/229 (%.1f%%) -> new %d/229 (%.1f%%)' % (
        (zk.p_replicate < 1e-3).sum(), 100 * (zk.p_replicate < 1e-3).mean(),
        (zk.q_excess < 0.05).sum(), 100 * (zk.q_excess < 0.05).mean()))
    print('  lowest penetrance_adj among newly-called regions: %.4f (old floor was 0.124)' % (
        cmp.loc[sig_new & ~sig_old, 'penetrance_adj'].min()))
    print('  median CI width on penetrance_adj: %.4f' % (cmp.penetrance_hi - cmp.penetrance_lo).median())

    summ = pd.DataFrame({
        'metric': ['n_donors_used', 'n_called', 'median_n_tested', 'zink_recovered',
                   'detection_floor_penetrance', 'median_ci_width'],
        'C07d_tier': [105, int(sig_old.sum()), float(cmp.n_tested_rep.median()),
                      int((zk.p_replicate < 1e-3).sum()), 0.124, np.nan],
        'C08b_adjusted': [int((b_i > 0.002).sum()), int(sig_new.sum()), float(cmp.n_tested_all.median()),
                          int((zk.q_excess < 0.05).sum()),
                          float(cmp.loc[sig_new, 'penetrance_adj'].min()),
                          float((cmp.penetrance_hi - cmp.penetrance_lo).median())]})
    summ.to_csv(OUT / 'penetrance_model_summary.tsv', sep='\t', index=False, float_format='%.5g')
    print('\nwrote', OUT / 'penetrance_adjusted.tsv.gz')


if __name__ == '__main__':
    main()
