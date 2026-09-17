#!/usr/bin/env python3
"""C06a_null_calibration.py -- empirical null for the read-level ASM test: split ONE haplotype's
reads at random into two halves and run the identical test. Any signal is false by construction,
so the resulting p-values give each donor's true calibration (and a per-donor inflation factor
that can be used to correct, or to exclude, donors).

Motivation (2026-09-17 chr20 spot check of the real calls): median genomic-inflation factor
lambda = 1.54 across 201 donors, mean 1.94, max 13.1, and lambda tracks donor PMD depth
(p < 1e-4) and XIST-skew clonality (r = 0.40) rather than chemistry or ancestry. Both are
mechanisms that make two random halves of one haplotype differ more than binomial sampling
implies: intermediate/regionally drifting methylation in PMDs, and less cell-to-cell averaging
in clonal lines. This script measures that directly.

Per (sample, chrom), for each haplotype separately: reads are shuffled and split 50/50, per-read
region fractions are computed exactly as in C02a, and ont_asm_caller.test_region_reads compares
the halves. Also reports, per region, the within-haplotype spread of read fractions, which is the
"dispersion" QC the user asked for.

Output: results/asm/null/<sample>_<chrom>.null.tsv.gz  (region_id, hap, n_reads_a, n_reads_b,
mu_a, mu_b, delta, pval, method, sd_reads) and <sample>_<chrom>.null_summary.tsv (per hap:
lambda_gc, frac p<0.05, mean |delta|, median read SD)

Usage: C06a_null_calibration.py <sample> <chrom> [seed]
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import comb

if not hasattr(math, 'comb'):      # allcools env is Python 3.7
    math.comb = lambda n, k: comb(n, k, exact=True)

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'asm'))
sys.path.insert(0, str(PROJDIR / 'github' / 'ont_asm_caller'))
from C02a_call_asm import DEFAULT_K, MIN_READS, REGIONS, hap_calls, read_fracs  # noqa: E402
from H02_filtered_locus_matrix import load_hmmflagger  # noqa: E402
from ont_asm_caller import test_region_reads  # noqa: E402

OUT = PROJDIR / 'results' / 'asm' / 'null'


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def main() -> None:
    sample, chrom = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    tag = f'{sample}_{chrom}'
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    hmm = load_hmmflagger(sample)
    reg = pd.read_csv(REGIONS, sep='\t')
    reg = reg[(reg.chrom == chrom) & (reg['set'] == 'cpg')].sort_values('start').reset_index(drop=True)

    rows, summ = [], []
    for hap in (1, 2):
        try:
            calls, _ = hap_calls(sample, hap, chrom, DEFAULT_K, hmm)
        except FileNotFoundError as err:
            log(f'# SKIP {tag} hap{hap}: {err}')
            continue
        if not len(calls):
            continue
        rids = calls.rid.unique()
        half = set(rng.choice(rids, size=len(rids) // 2, replace=False).tolist())
        side = calls.rid.isin(half)
        fr_a = read_fracs(calls[side], reg)
        fr_b = read_fracs(calls[~side], reg)
        by_a = {i: g.frac.values for i, g in fr_a.groupby('reg')}
        by_b = {i: g.frac.values for i, g in fr_b.groupby('reg')}
        n_hap = 0
        for i in sorted(set(by_a) & set(by_b)):
            a, b = by_a[i], by_b[i]
            if len(a) < MIN_READS or len(b) < MIN_READS:
                continue
            res = test_region_reads(chrom, int(reg.start[i]), int(reg.end[i]), a[:, None], b[:, None],
                                    min_calls=1, min_reads=MIN_READS)
            if res is None:
                continue
            rows.append((reg.region_id[i], hap, res.n_reads1, res.n_reads2, res.mu1_hat, res.mu2_hat,
                         res.delta, res.pval, res.method, float(np.std(np.r_[a, b], ddof=1))))
            n_hap += 1
        p = np.array([r[7] for r in rows[-n_hap:]]) if n_hap else np.array([])
        if len(p):
            chi2 = stats.chi2.isf(np.clip(p, 1e-300, 1), 1)
            summ.append({'sample': sample, 'chrom': chrom, 'hap': hap, 'n_regions': n_hap,
                         'lambda_null': float(np.median(chi2) / 0.4549),
                         'frac_p05': float((p < 0.05).mean()), 'frac_p001': float((p < 0.001).mean()),
                         'mean_absdelta': float(np.mean([abs(r[6]) for r in rows[-n_hap:]])),
                         'median_read_sd': float(np.median([r[9] for r in rows[-n_hap:]]))})
            log(f'{tag} hap{hap}: {n_hap} regions, lambda_null '
                f'{summ[-1]["lambda_null"]:.2f}, frac p<0.05 {summ[-1]["frac_p05"]:.3f}')

    cols = ['region_id', 'hap', 'n_reads_a', 'n_reads_b', 'mu_a', 'mu_b', 'delta', 'pval', 'method', 'sd_reads']
    pd.DataFrame(rows, columns=cols).to_csv(OUT / f'{tag}.null.tsv.gz', sep='\t', index=False,
                                            float_format='%.6g')
    pd.DataFrame(summ).to_csv(OUT / f'{tag}.null_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
