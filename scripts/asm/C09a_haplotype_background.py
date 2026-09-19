#!/usr/bin/env python3
"""C09a_haplotype_background.py -- does ASM penetrance depend on WHICH local haplotype carries the
lead allele, rather than on the lead allele alone?

The question (RESULTS §13). Among donors heterozygous at a `genotype_linked` locus's lead variant
-- i.e. carriers of the putative causal allele -- median penetrance is only 0.525, and just 2.1%
of loci are fully penetrant. That incomplete penetrance is locus-intrinsic (locus identity 18.1%
of the variance vs donor identity 3.8%) and is not read depth. The leading explanation is that
the lead SNV is a TAG, not the cause: it sits on more than one local haplotype, and only some of
those haplotypes carry whatever actually drives the ASM. If so, penetrance should be predicted by
the haplotype background in phase with the lead allele.

This matters because it is the whole point of an ancestrally diverse cohort: haplotype
frequencies differ strongly between superpopulations even where tag-allele frequencies do not, so
a haplotype-background mechanism would explain the ancestry-differential penetrance signal
(7.6% of loci at p<0.05 vs 5% expected) without invoking anything trans.

Design, per locus, restricted to donors heterozygous at the lead variant:
  - identify which haplotype (hap1 or hap2) carries the lead ALT allele in each such donor;
  - read every common background variant (AF >= MIN_AF) within +-WINDOW bp *in phase with that
    allele* -- i.e. the allele on the ALT-carrying haplotype. This is the key step: it converts a
    diploid genotype into the single haplotype that is doing the work;
  - two tests of ASM status against that background:
      `best`     per-variant scan, Fisher/hypergeometric, minimum p over background variants;
      `cluster`  split carriers in two on PC1 of the phased background matrix, then one test.
  - **both are calibrated by permutation from the start** (shuffle ASM among lead-het carriers,
    recompute the statistic, B times). The C08c lesson: a naive Bonferroni over correlated
    variants was badly wrong in both directions, so no threshold here is trusted analytically.

A locus where `p_perm_best` is small but the lead variant alone already explained ASM is not
interesting; the informative quantity is the lift, so the lead-only penetrance split is carried
alongside for comparison.

NOTE -- the lead variant is RE-DERIVED here, not read from C07c/C07d. Both of those write their
tables with `float_format='%.5g'`, which truncates `lead_pos` to five significant figures
(606320 is stored as `6.0632e+05`), so the recorded positions do not identify a variant. The
analyses in those scripts are unaffected because they held the position in memory; only the
on-disk column is corrupt. Re-deriving also gives a free cross-check: `lead_p_rederived` here
should reproduce the stored `lead_p`.

Output: results/asm/model/hap_background/<chrom>.hapbg.tsv.gz
  region_id, n_het, n_het_asm, pen_het, n_bg_variants, best_p, best_pos, best_dist,
  best_pen_hi/lo (penetrance in the two background groups), p_perm_best, cluster_p,
  cluster_pen_hi/lo, p_perm_cluster, pc1_var

Usage: C09a_haplotype_background.py <chrom> [n_perm]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import hypergeom as _hg

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'asm'))
REP = PROJDIR / 'results' / 'asm' / 'replication'
OUT = PROJDIR / 'results' / 'asm' / 'model' / 'hap_background'
WINDOW = 25_000
FLANK = 5_000          # lead-variant search radius, matches C07c
MIN_AF = 0.05
MIN_HET = 20
MIN_GRP = 5
CALIBRATED = ('discovery', 'replication')


def pval_table(N: int, K: int, nvals: np.ndarray) -> np.ndarray:
    """Two-sided hypergeometric p for every (n, a): 2 * min(upper tail, lower tail), clipped.

    Computed directly from log-binomials rather than scipy.stats.hypergeom, whose broadcasting
    path is unreliable in the allcools env's scipy (it drops into a scalar branch and raises
    "truth value of an array is ambiguous" on some shapes). Fully vectorised over (n, a).
    """
    a = np.arange(int(K) + 1)                       # successes drawn
    nv = np.asarray(nvals, dtype=np.int64)          # draws
    lg = gammaln(np.arange(N + 2))                  # lg[i] = log((i-1)!)

    def lchoose(n, k):
        n = np.asarray(n)
        k = np.asarray(k)
        bad = (k < 0) | (k > n) | (n < 0)
        kk = np.clip(k, 0, np.maximum(n, 0))
        out = lg[n + 1] - lg[kk + 1] - lg[np.maximum(n - kk, 0) + 1]
        return np.where(bad, -np.inf, out)

    A = a[None, :]
    Nv = nv[:, None]
    logpmf = lchoose(K, A) + lchoose(N - K, Nv - A) - lchoose(N, Nv)
    pmf = np.exp(logpmf)
    pmf = np.where(np.isfinite(logpmf), pmf, 0.0)
    upper = np.cumsum(pmf[:, ::-1], axis=1)[:, ::-1]   # P(X >= a)
    lower = np.cumsum(pmf, axis=1)                     # P(X <= a)
    return np.clip(2 * np.minimum(upper, lower), 0, 1)


def main() -> None:
    chrom = sys.argv[1]
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / ('%s.hapbg.tsv.gz' % chrom)
    if dst.exists():
        print('exists, skip:', dst)
        return

    from C07c_candidate_genotypes import load_genotypes, VCF, BCFTOOLS  # noqa: E402
    import subprocess

    tiers = pd.read_csv(REP / 'donor_tiers.tsv', sep='\t')
    cand = pd.read_csv(REP / 'candidates_replication.tsv.gz', sep='\t',
                       usecols=['region_id', 'chrom', 'start', 'end', 'rep_class', 'lead_p'],
                       low_memory=False)
    cand = cand[(cand.chrom == chrom) & (cand.rep_class == 'genotype_linked')]
    cand = cand.sort_values('start').reset_index(drop=True)
    print('%s: %d genotype_linked loci with a lead variant, B=%d' % (chrom, len(cand), B), flush=True)
    if cand.empty:
        pd.DataFrame().to_csv(dst, sep='\t', index=False)
        return

    vcf_samples = subprocess.run([BCFTOOLS, 'query', '-l', str(VCF)], capture_output=True,
                                 text=True, check=True).stdout.split()
    sidx = {s: i for i, s in enumerate(vcf_samples)}

    w = np.stack([np.maximum(cand.start.values - WINDOW, 0), cand.end.values + WINDOW], axis=1)
    merged = [list(w[0])]
    for s, e in w[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    pos1, ref, alt, h1, h2 = load_genotypes(chrom, np.array(merged, dtype=np.int64), vcf_samples)
    p0 = pos1 - 1
    af = (h1.sum(1) + h2.sum(1)) / (2.0 * h1.shape[1]) if len(pos1) else np.zeros(0)
    print('%s: %d SNVs loaded' % (chrom, len(pos1)), flush=True)

    cset = set(cand.region_id)
    calls = []
    for s in tiers['sample']:
        f = REP / 'calls' / ('%s.tsv.gz' % s)
        if f.exists():
            d = pd.read_csv(f, sep='\t', usecols=['sample', 'region_id', 'asm_rep'])
            calls.append(d[d.region_id.isin(cset)])
    calls = pd.concat(calls, ignore_index=True).merge(tiers[['sample', 'tier']], on='sample')
    calls = calls[calls['sample'].isin(sidx) & calls.tier.isin(CALIBRATED)]
    by_region = {r: g for r, g in calls.groupby('region_id')}

    rng = np.random.default_rng(1)
    rows = []
    drop = {k: 0 for k in ('no_calls','lead_missing','few_het','few_grp','few_bg','few_ok','svd')}
    for ii, r in enumerate(cand.itertuples(index=False)):
        g = by_region.get(r.region_id)
        if g is None:
            drop['no_calls'] += 1
            continue
        cols = np.array([sidx[s] for s in g['sample']])
        asm_all = g.asm_rep.values.astype(bool)

        # --- re-derive the lead variant exactly as C07c does (its stored lead_pos is corrupt) ---
        lo5, hi5 = np.searchsorted(p0, [r.start - FLANK, r.end + FLANK])
        v5 = np.arange(lo5, hi5)
        if len(v5) == 0:
            drop['lead_missing'] += 1
            continue
        Hall = (h1[np.ix_(v5, cols)] != h2[np.ix_(v5, cols)])
        Nall, Kall = len(cols), int(asm_all.sum())
        nn = Hall.sum(1)
        okl = (nn >= 2) & (Nall - nn >= 2)
        if not okl.any() or Kall < 1:
            drop['lead_missing'] += 1
            continue
        aa = (Hall & asm_all[None, :]).sum(1)
        pl = np.where(okl, _hg.sf(aa - 1, Nall, Kall, nn), 1.0)
        jl = int(np.argmin(pl))
        lead_i = int(v5[jl])
        lead_p_red = float(pl[jl])

        l1, l2 = h1[lead_i, cols], h2[lead_i, cols]
        het = l1 != l2
        if het.sum() < MIN_HET:
            drop['few_het'] += 1
            continue
        asm = asm_all[het]
        N, K = int(het.sum()), int(asm.sum())
        if K < MIN_GRP or N - K < MIN_GRP:
            drop['few_grp'] += 1
            continue
        hc = cols[het]
        alt_on_h1 = l1[het] == 1

        lo, hi = np.searchsorted(p0, [r.start - WINDOW, r.end + WINDOW])
        sel = np.arange(lo, hi)
        sel = sel[(af[sel] >= MIN_AF) & (af[sel] <= 1 - MIN_AF) & (sel != lead_i)]
        if len(sel) < 2:
            drop['few_bg'] += 1
            continue
        # allele carried ON the lead-ALT haplotype  (V x N_het)
        P = np.where(alt_on_h1[None, :], h1[np.ix_(sel, hc)], h2[np.ix_(sel, hc)]).astype(np.int8)
        n = P.sum(1)
        ok = (n >= MIN_GRP) & (N - n >= MIN_GRP)
        if not ok.any():
            drop['few_ok'] += 1
            continue
        P, sel, n = P[ok], sel[ok], n[ok]
        nun, ninv = np.unique(n, return_inverse=True)
        table = pval_table(N, K, nun)

        a_obs = (P & asm[None, :]).sum(1)
        p_obs = table[ninv, a_obs]
        j = int(np.argmin(p_obs))
        best_p = float(p_obs[j])
        grp = P[j].astype(bool)
        pen_hi = float(asm[grp].mean())
        pen_lo = float(asm[~grp].mean())
        if pen_lo > pen_hi:
            pen_hi, pen_lo = pen_lo, pen_hi

        # PC1 cluster of the phased background
        X = P.astype(np.float64).T
        X -= X.mean(0)
        try:
            pc1 = np.linalg.svd(X, full_matrices=False)[0][:, 0]
            sv = np.linalg.svd(X, full_matrices=False)[1]
            pc1_var = float(sv[0] ** 2 / max((sv ** 2).sum(), 1e-9))
        except np.linalg.LinAlgError:
            drop['svd'] += 1
            continue
        cl = pc1 > np.median(pc1)
        nc = int(cl.sum())
        cl_p, cl_hi, cl_lo = np.nan, np.nan, np.nan
        if MIN_GRP <= nc <= N - MIN_GRP:
            tb = pval_table(N, K, np.array([nc]))
            cl_p = float(tb[0, int((cl & asm).sum())])
            cl_hi, cl_lo = float(asm[cl].mean()), float(asm[~cl].mean())
            if cl_lo > cl_hi:
                cl_hi, cl_lo = cl_lo, cl_hi

        # permutations: shuffle which carriers have ASM
        Ap = np.zeros((N, B), dtype=np.int8)
        for b in range(B):
            Ap[rng.choice(N, K, replace=False), b] = 1
        cnt = (P.astype(np.float32) @ Ap.astype(np.float32)).astype(np.int64)
        perm_best = table[ninv[:, None], cnt].min(0)
        p_perm_best = (1 + int((perm_best <= best_p).sum())) / (B + 1.0)
        p_perm_cl = np.nan
        if np.isfinite(cl_p):
            cc = (cl.astype(np.float32)[None, :] @ Ap.astype(np.float32))[0].astype(np.int64)
            tb = pval_table(N, K, np.array([nc]))[0]
            p_perm_cl = (1 + int((tb[cc] <= cl_p).sum())) / (B + 1.0)

        rows.append({'region_id': r.region_id, 'chrom': chrom, 'lead_pos': int(pos1[lead_i]),
                     'lead_p_rederived': lead_p_red, 'lead_p_stored': float(r.lead_p),
                     'n_het': N, 'n_het_asm': K, 'pen_het': K / N, 'n_bg_variants': int(len(sel)),
                     'best_p': best_p, 'best_pos': int(pos1[sel[j]]),
                     'best_dist': int(abs(int(pos1[sel[j]]) - int(pos1[lead_i]))),
                     'best_pen_hi': pen_hi, 'best_pen_lo': pen_lo, 'p_perm_best': p_perm_best,
                     'cluster_p': cl_p, 'cluster_pen_hi': cl_hi, 'cluster_pen_lo': cl_lo,
                     'p_perm_cluster': p_perm_cl, 'pc1_var': pc1_var})
        if (ii + 1) % 250 == 0:
            print('  %d/%d' % (ii + 1, len(cand)), flush=True)

    print('drop reasons:', drop, flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(dst, sep='\t', index=False, float_format='%.5g')
    if len(out):
        print('\n%s: %d loci tested' % (chrom, len(out)))
        print('  median pen_het %.3f ; median background variants %d' % (
            out.pen_het.median(), out.n_bg_variants.median()))
        for col in ('p_perm_best', 'p_perm_cluster'):
            v = out[col].dropna()
            if len(v):
                print('  %s: frac<0.05 %.3f (null 0.05), frac<0.01 %.3f (null 0.01)' % (
                    col, (v < 0.05).mean(), (v < 0.01).mean()))
        sig = out[out.p_perm_best < 0.05]
        if len(sig):
            print('  among permutation-significant loci: penetrance %.3f vs %.3f between'
                  ' background groups (lift %.3f)' % (
                      sig.best_pen_hi.median(), sig.best_pen_lo.median(),
                      (sig.best_pen_hi - sig.best_pen_lo).median()))
    print('wrote', dst)


if __name__ == '__main__':
    main()
