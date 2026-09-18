#!/usr/bin/env python3
"""C07b_retest_candidates.py -- ASM replication, step 2: re-test the C07a candidates in every
donor, with per-donor genomic control.

For one donor: read its C02a calls on all autosomes, keep the candidate regions, then
  - genomic control: chi2_1df(p) / lambda_gc -> p_gc (lambda from C07a donor_tiers.tsv; only
    applied when lambda > 1 -- deflated donors are left as-is, the test is already conservative)
  - BH across the candidate set only (q_cand): candidates were chosen in *other* donors, so the
    multiple-testing burden in a replication donor is the candidate count, not the genome
  - asm_rep: q_cand < Q_MAX and |delta| >= MIN_DELTA (same effect-size rule as C04a)
  - asm_gw: the donor's original genome-wide call (C04a), for comparison
Candidates the donor never tested (too few reads per haplotype) are simply absent.

Discovery-tier donors are re-tested the same way, but their own calls chose the candidates, so
C07d computes penetrance on replication-tier donors (and reports discovery separately).

Output: results/asm/replication/calls/<sample>.tsv.gz
  sample, region_id, n_reads1, n_reads2, delta (hap1 - hap2), pval, p_gc, q_cand, asm_rep, asm_gw

Usage: C07b_retest_candidates.py <sample>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
CALLS = PROJDIR / 'results' / 'asm' / 'calls'
SIG = PROJDIR / 'results' / 'asm' / 'genome' / 'per_donor'
REP = PROJDIR / 'results' / 'asm' / 'replication'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
Q_MAX = 0.05
MIN_DELTA = 0.2


def main(sample: str) -> None:
    out = REP / 'calls' / f'{sample}.tsv.gz'
    out.parent.mkdir(parents=True, exist_ok=True)
    tiers = pd.read_csv(REP / 'donor_tiers.tsv', sep='\t').set_index('sample')
    lam = tiers.lambda_gc.get(sample, np.nan)
    cand = set(pd.read_csv(REP / 'candidates.tsv.gz', sep='\t', usecols=['region_id']).region_id)

    cols = ['region_id', 'n_reads1', 'n_reads2', 'delta', 'pval']
    parts = []
    for c in AUTOSOMES:
        f = CALLS / f'{sample}_{c}.asm.tsv.gz'
        if not f.exists():
            print(f'WARNING: missing {f.name}', flush=True)
            continue
        d = pd.read_csv(f, sep='\t', usecols=cols)
        parts.append(d[d.region_id.isin(cand)])
    d = pd.concat(parts, ignore_index=True)

    p = d.pval.clip(lower=1e-300).values
    if np.isfinite(lam) and lam > 1:
        d['p_gc'] = chi2.sf(chi2.isf(p, 1) / lam, 1)
    else:
        d['p_gc'] = p
    d['q_cand'] = multipletests(d.p_gc.values, method='fdr_bh')[1]
    d['asm_rep'] = (d.q_cand < Q_MAX) & (d.delta.abs() >= MIN_DELTA)
    sig = pd.read_csv(SIG / f'{sample}.asm_sig.tsv.gz', sep='\t', usecols=['region_id'])
    d['asm_gw'] = d.region_id.isin(set(sig.region_id))
    d.insert(0, 'sample', sample)
    d.to_csv(out, sep='\t', index=False, float_format='%.6g')
    print(f'{sample} (lambda_gc {lam:.2f}, tier {tiers.tier.get(sample, "NA")}): '
          f'{len(d):,}/{len(cand):,} candidates tested, {int(d.asm_rep.sum()):,} asm_rep, '
          f'{int(d.asm_gw.sum()):,} asm_gw', flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
