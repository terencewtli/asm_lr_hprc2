#!/usr/bin/env python3
"""C08c_lead_permutation.py -- permutation null for the lead-variant scan, per chromosome.

Why. C07c picks each region's lead variant by `argmin` of a hypergeometric p over every het SNV
within +-5 kb (median 75, mean 83, max 1,219 variants), and C07d then thresholds that minimum at
an UNCORRECTED 1e-4. The reported lead_p is therefore a minimum-of-many-tests statistic being
read as a single-test p. Naive upper bound: ~986 spurious leads across 118,796 candidates, i.e.
~5% of the 18,225 regions passing the lead filter. LD makes the tests dependent, so the real
number needs a permutation, not a Bonferroni guess.

Null. Hold the variant set, the het matrix and K = number of ASM-calling donors fixed; permute
WHICH calibrated donors carry the ASM call. That destroys any variant-ASM association while
preserving region-level penetrance, per-variant allele frequency and the full LD structure --
so the permuted min-p distribution is the correct null for "best of these variants".

Implementation. hypergeom.sf(a-1, N, K, n) depends only on (a, n) once N and K are fixed, so the
p-values are read from a precomputed (a, n) lookup and each region's B permutations reduce to one
BLAS call: Hc (V x N) @ Aperm (N x B) -> counts (V x B).

Empirical p for the region = (1 + #{min_p_perm <= min_p_obs}) / (B + 1).

Output: results/asm/model/lead_perm/<chrom>.lead_perm.tsv.gz
  region_id, lead_p (as C07c), n_snv_tested_5kb, p_perm, and the permuted-min-p quantiles.

Usage: C08c_lead_permutation.py <chrom> [n_perm]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'asm'))
REP = PROJDIR / 'results' / 'asm' / 'replication'
OUT = PROJDIR / 'results' / 'asm' / 'model' / 'lead_perm'
FLANK = 5000
NEAR_MAX = 50_000
CALIBRATED = ('discovery', 'replication')


def main() -> None:
    chrom = sys.argv[1]
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / ('%s.lead_perm.tsv.gz' % chrom)
    if dst.exists():
        print('exists, skip:', dst)
        return

    from C07c_candidate_genotypes import load_genotypes, VCF, BCFTOOLS  # noqa: E402
    import subprocess

    tiers = pd.read_csv(REP / 'donor_tiers.tsv', sep='\t')
    cand = pd.read_csv(REP / 'candidates.tsv.gz', sep='\t',
                       usecols=['region_id', 'set', 'chrom', 'start', 'end'])
    cand = cand[cand.chrom == chrom].sort_values('start').reset_index(drop=True)
    print('%s: %d candidates, B=%d' % (chrom, len(cand), B), flush=True)

    vcf_samples = subprocess.run([BCFTOOLS, 'query', '-l', str(VCF)], capture_output=True,
                                 text=True, check=True).stdout.split()
    sidx = {s: i for i, s in enumerate(vcf_samples)}

    w = np.stack([np.maximum(cand.start.values - NEAR_MAX, 0), cand.end.values + NEAR_MAX], axis=1)
    merged = [list(w[0])]
    for s, e in w[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    pos1, ref, alt, h1, h2 = load_genotypes(chrom, np.array(merged, dtype=np.int64), vcf_samples)
    p0 = pos1 - 1
    het = h1 != h2
    print('%s: %d SNVs loaded' % (chrom, len(pos1)), flush=True)

    cset = set(cand.region_id)
    calls = []
    for s in tiers['sample']:
        f = REP / 'calls' / ('%s.tsv.gz' % s)
        if f.exists():
            d = pd.read_csv(f, sep='\t', usecols=['sample', 'region_id', 'asm_rep'])
            calls.append(d[d.region_id.isin(cset)])
    calls = pd.concat(calls, ignore_index=True).merge(tiers[['sample', 'tier']], on='sample')
    calls = calls[calls['sample'].isin(sidx)]
    by_region = {r: g for r, g in calls.groupby('region_id')}

    rng = np.random.default_rng(1)
    rows = []
    for ii, r in enumerate(cand.itertuples(index=False)):
        g = by_region.get(r.region_id)
        if g is None:
            continue
        lo, hi = np.searchsorted(p0, [r.start - NEAR_MAX, r.end + NEAR_MAX])
        vp = p0[lo:hi]
        dist = np.maximum(0, np.maximum(r.start - vp, vp - (r.end - 1)))
        in5 = dist <= FLANK
        cols = np.array([sidx[s] for s in g['sample']])
        cal = g.tier.isin(CALIBRATED).values
        asm = g.asm_rep.values.astype(bool)
        N, K = int(cal.sum()), int((asm & cal).sum())
        if N < 4 or K < 1 or not in5.any():
            continue
        H = het[lo:hi][:, cols]
        Hc = H[np.flatnonzero(in5)][:, cal]
        n = Hc.sum(1)
        ok = (n >= 2) & (N - n >= 2)
        if not ok.any():
            continue
        Hc = Hc[ok]
        n = n[ok]
        V = Hc.shape[0]

        # p lookup over (a, n): a in 0..K, n in observed set
        nun, ninv = np.unique(n, return_inverse=True)
        av = np.arange(0, K + 1)
        table = hypergeom.sf(av[None, :] - 1, N, K, nun[:, None])       # (len(nun), K+1)

        a_obs = (Hc & asm[cal][None, :]).sum(1)
        p_obs = table[ninv, a_obs]
        min_obs = float(p_obs.min())

        # B permutations of which calibrated donors call ASM, K preserved
        Ap = np.zeros((N, B), dtype=np.float32)
        for b in range(B):
            Ap[rng.choice(N, K, replace=False), b] = 1.0
        counts = (Hc.astype(np.float32) @ Ap).astype(np.int64)           # (V, B)
        p_perm = table[ninv[:, None], counts]
        min_perm = p_perm.min(0)

        rows.append({'region_id': r.region_id, 'set': r.set, 'n_snv_tested_5kb': int(V),
                     'N_cal': N, 'K_cal': K, 'lead_p': min_obs,
                     'p_perm': (1 + int((min_perm <= min_obs).sum())) / (B + 1.0),
                     'perm_min_p05': float(np.quantile(min_perm, 0.05)),
                     'perm_min_median': float(np.median(min_perm))})
        if (ii + 1) % 1000 == 0:
            print('  %d/%d' % (ii + 1, len(cand)), flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(dst, sep='\t', index=False, float_format='%.5g')
    if len(out):
        hit = out.lead_p < 1e-4
        print('\n%s: %d regions scanned; %d with lead_p<1e-4' % (chrom, len(out), hit.sum()))
        if hit.any():
            print('  of those, permutation p >= 0.05 (i.e. NOT significant once corrected): %d (%.1f%%)'
                  % ((out.loc[hit, 'p_perm'] >= 0.05).sum(), 100 * (out.loc[hit, 'p_perm'] >= 0.05).mean()))
        print('  median permuted min-p: %.3g (this is the real per-region threshold scale)'
              % out.perm_min_median.median())
    print('wrote', dst)


if __name__ == '__main__':
    main()
