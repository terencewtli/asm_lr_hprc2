#!/usr/bin/env python3
"""C10a_epiallele_nme.py -- epiallele entropy (CPEL's NME) per donor at the replicating ASM loci.

SCOPE, deliberately narrow (see JOURNAL 2026-09-19). This is NOT a second ASM discovery axis --
no loci are called here. It answers one pre-registered question that RESULTS §13 cannot:

    does epiallele entropy predict penetrance?

§13 establishes that incomplete penetrance is locus-intrinsic (locus identity 18.1% of the
variance vs donor 3.8%, depth-independent) but every structural predictor of `pen_het` came back
at |rho| <= 0.235. NME is the one candidate the project has not measured. The distinction it
draws: is a low-penetrance locus a SWITCH (some donors show clean two-state ASM, others none) or
a GRADED element (reads individually intermediate in everyone)? The read-level caller is blind to
this by construction -- `readlevel.test_region_reads` collapses each molecule to a fraction, so a
haplotype pair with equal means and different disorder is invisible to it.

There is also a directional prediction to test rather than a fishing expedition: `asm_lr`'s W12
pilot (one donor) found genotype-proximal ASM is MORE deterministic -- lower within-haplotype
entropy, p < 1e-54. Either that replicates at n=174 or it does not.

What is computed, per (donor, locus):
  nme1, nme2   normalised methylation entropy (bits/CpG) of each haplotype's fitted CPEL chain
  t_nme        |nme1 - nme2|, the entropy-imbalance statistic
  mml1, mml2   mean methylation level of each chain -- free, and a sanity check against the
               caller's own delta (they should agree closely; if not, the fit is wrong)
  nme_null     the SAME statistic under the empirical null: hap1's reads split in half by WHOLE
               MOLECULE and the two halves compared. Built in from the start, not bolted on --
               `pattern.py` has a simulation benchmark but NO real-data validation, so no
               locus-level claim is possible without it.

NOT computed: T_PDM. `jsd_chains` is Monte-Carlo over the fitted distributions, it is the
expensive statistic, and it needs the stratified-null machinery in `null.py` to be interpretable.
MML/NME answer the §13 question; PDM is a different paper.

Outputs: results/asm/model/epiallele/<sample>.nme.tsv.gz

Usage: C10a_epiallele_nme.py <sample> [--limit N]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'asm'))
sys.path.insert(0, str(PROJDIR / 'github' / 'ont_asm_caller'))
REP = PROJDIR / 'results' / 'asm' / 'replication'
OUT = PROJDIR / 'results' / 'asm' / 'model' / 'epiallele'
AUTOSOMES = ['chr%d' % i for i in range(1, 23)]
MIN_READS = 8          # a 2-state fit on fewer molecules is fitting noise (cf. H04's MIN_READS)
MIN_CPG = 5
REPLICATING = ('imprinting', 'genotype_linked', 'genotype_indep', 'other')


def build_matrix(calls: pd.DataFrame, pos_index: dict) -> np.ndarray:
    """read x CpG matrix of 0/1/nan from long-form (rid, pos, call) rows of one region."""
    rids = calls.rid.values
    uniq, rinv = np.unique(rids, return_inverse=True)
    m = np.full((len(uniq), len(pos_index)), np.nan, dtype=np.float64)
    cidx = calls.pos.map(pos_index).values
    ok = ~pd.isna(cidx)
    m[rinv[ok], cidx[ok].astype(int)] = calls.call.values[ok]
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('sample')
    ap.add_argument('--limit', type=int, default=0, help='stop after N loci (benchmarking)')
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    dst = OUT / ('%s.nme.tsv.gz' % a.sample)
    if dst.exists():
        print('exists, skip:', dst)
        return

    from C02a_call_asm import DEFAULT_K, REGIONS, hap_calls  # noqa: E402
    from H02_filtered_locus_matrix import load_hmmflagger    # noqa: E402
    from ont_asm_caller import fit_chain, chain_entropy_bits, chain_mml, split_reads  # noqa: E402

    cand = pd.read_csv(REP / 'candidates_replication.tsv.gz', sep='\t',
                       usecols=['region_id', 'chrom', 'start', 'end', 'rep_class'], low_memory=False)
    cand = cand[cand.rep_class.isin(REPLICATING)]
    print('%s: %d replicating loci' % (a.sample, len(cand)), flush=True)

    reg_all = pd.read_csv(REGIONS, sep='\t')
    hmm = load_hmmflagger(a.sample)
    rng = np.random.default_rng(1)

    def nme_of(m):
        if np.isfinite(m).any(axis=1).sum() < MIN_READS:
            return None
        h, J, _, _ = fit_chain(m)
        return chain_entropy_bits(h, J) / m.shape[1], chain_mml(h, J)

    rows, done = [], 0
    for chrom in AUTOSOMES:
        sub = cand[cand.chrom == chrom]
        if sub.empty:
            continue
        try:
            c1, _ = hap_calls(a.sample, 1, chrom, DEFAULT_K, hmm)
            c2, _ = hap_calls(a.sample, 2, chrom, DEFAULT_K, hmm)
        except Exception as e:                      # missing chain / modbed for this chrom
            print('  %s: skip (%s: %s)' % (chrom, type(e).__name__, e), flush=True)
            continue
        if c1.empty or c2.empty:
            continue
        regc = reg_all[(reg_all.chrom == chrom) & (reg_all['set'] == 'cpg')]
        starts = regc.set_index('region_id')[['start', 'end']]
        c1 = c1.sort_values('pos'); c2 = c2.sort_values('pos')
        p1, p2 = c1.pos.values, c2.pos.values

        for r in sub.itertuples(index=False):
            if r.region_id not in starts.index:
                continue
            lo1, hi1 = np.searchsorted(p1, [r.start, r.end])
            lo2, hi2 = np.searchsorted(p2, [r.start, r.end])
            g1, g2 = c1.iloc[lo1:hi1], c2.iloc[lo2:hi2]
            if len(g1) == 0 or len(g2) == 0:
                continue
            pos_u = np.union1d(g1.pos.unique(), g2.pos.unique())
            if len(pos_u) < MIN_CPG:
                continue
            pidx = {p: i for i, p in enumerate(pos_u)}
            m1, m2 = build_matrix(g1, pidx), build_matrix(g2, pidx)
            s1, s2 = nme_of(m1), nme_of(m2)
            if s1 is None or s2 is None:
                continue
            # empirical null: split hap1's reads by WHOLE MOLECULE and compare the halves
            nme_null = np.nan
            if m1.shape[0] >= 2 * MIN_READS:
                ma, mb = split_reads(rng, m1, m1)
                sa, sb = nme_of(ma), nme_of(mb)
                if sa is not None and sb is not None:
                    nme_null = abs(sa[0] - sb[0])
            rows.append({'sample': a.sample, 'region_id': r.region_id, 'rep_class': r.rep_class,
                         'n_cpg': len(pos_u), 'n_reads1': int(np.isfinite(m1).any(1).sum()),
                         'n_reads2': int(np.isfinite(m2).any(1).sum()),
                         'nme1': s1[0], 'nme2': s2[0], 't_nme': abs(s1[0] - s2[0]),
                         'mml1': s1[1], 'mml2': s2[1], 't_mml': abs(s1[1] - s2[1]),
                         'nme_null': nme_null})
            done += 1
            if a.limit and done >= a.limit:
                break
        print('  %s done (%d loci so far)' % (chrom, done), flush=True)
        if a.limit and done >= a.limit:
            break

    out = pd.DataFrame(rows)
    out.to_csv(dst, sep='\t', index=False, float_format='%.5g')
    if len(out):
        print('\n%s: %d loci' % (a.sample, len(out)))
        print('  NME per hap: median %.3f (IQR %.3f-%.3f)' % (
            np.median(np.r_[out.nme1, out.nme2]),
            np.quantile(np.r_[out.nme1, out.nme2], .25), np.quantile(np.r_[out.nme1, out.nme2], .75)))
        print('  t_nme  observed median %.4f  vs  null median %.4f' % (
            out.t_nme.median(), out.nme_null.median()))
        print('  t_mml  median %.4f   (sanity: should track the caller delta)' % out.t_mml.median())
        print('  by class:\n' + out.groupby('rep_class')[['t_nme', 'nme_null']].median().round(4).to_string())
    print('wrote', dst)


if __name__ == '__main__':
    main()
