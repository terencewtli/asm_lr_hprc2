#!/usr/bin/env python3
"""A03b_chr20_pmd_survey_summary.py -- aggregate A03a's per-(sample, hap) chr20 outputs.

Reference PMD set = hg38 chr20 CpGs inside a PMD on BOTH NA19338 haplotypes; reference
non-PMD = CpGs outside PMDs on both (discordant CpGs dropped). For every donor/hap:
  - mean methylation inside vs outside the reference PMDs (the "is it a caller limitation?"
    test: if even high-global-methylation donors are lower inside, the domains exist there too,
    just too shallow for dnmtools pmd to call)
  - its own PMD burden (from its own dnmtools pmd call)
Per donor, hap1 vs hap2:
  - mean |hap1 - hap2| over 10kb hg38 bins (>=10 shared CpGs), bins split by reference label
  - own-PMD call concordance between haps (Jaccard over shared CpGs)

Outputs (results/qc/data/chr20_pmd_survey/): per_hap.tsv, per_donor.tsv, by_superpop.tsv

Usage: A03b_chr20_pmd_survey_summary.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
SURVEY = PROJDIR / 'results' / 'qc' / 'data' / 'chr20_pmd_survey'
PER_HAP = SURVEY / 'per_hap'
MANIFEST = PROJDIR / 'tsv' / 'meta' / 'hprc2_sample_manifest.tsv'
GLOBAL_METH = PROJDIR / 'results' / 'qc' / 'data' / 'global_methylation_covariates_full_cohort.tsv'
REF = 'NA19338'
BIN = 10_000
MIN_BIN_CPG = 10


def load_cpg(sample: str, hap: int) -> pd.DataFrame:
    return pd.read_csv(PER_HAP / f'{sample}_hap{hap}.hg38_cpg.tsv.gz', sep='\t')


def main() -> None:
    summ = pd.concat([pd.read_csv(f, sep='\t') for f in sorted(PER_HAP.glob('*.summary.tsv'))],
                     ignore_index=True)
    print(f'{len(summ)} (sample, hap) summaries, {summ["sample"].nunique()} donors', flush=True)

    r1, r2 = load_cpg(REF, 1), load_cpg(REF, 2)
    ref = r1.merge(r2, on='pos', suffixes=('_1', '_2'))
    ref = ref[ref.in_own_pmd_1 == ref.in_own_pmd_2][['pos', 'in_own_pmd_1']]
    ref = ref.rename(columns={'in_own_pmd_1': 'ref_pmd'}).set_index('pos').ref_pmd
    print(f'reference CpGs: {len(ref):,} ({ref.mean():.3f} in PMD on both {REF} haps)', flush=True)

    rows = []
    donor_rows = []
    for sample, grp in summ.groupby('sample'):
        tabs = {}
        for hap in grp.hap:
            t = load_cpg(sample, hap)
            t['ref_pmd'] = t.pos.map(ref)
            t = t.dropna(subset=['ref_pmd'])
            tabs[hap] = t
            rows.append({'sample': sample, 'hap': hap,
                         'meth_in_ref_pmd': t.frac[t.ref_pmd == 1].mean(),
                         'meth_out_ref_pmd': t.frac[t.ref_pmd == 0].mean(),
                         'n_cpg_ref': len(t)})
        if len(tabs) < 2:
            continue
        m = tabs[1].merge(tabs[2], on='pos', suffixes=('_1', '_2'))
        m['bin'] = m.pos // BIN
        b = m.groupby('bin').agg(n=('pos', 'size'), f1=('frac_1', 'mean'), f2=('frac_2', 'mean'),
                                 ref_pmd=('ref_pmd_1', 'mean'))
        b = b[b.n >= MIN_BIN_CPG]
        b['absdiff'] = (b.f1 - b.f2).abs()
        both = ((m.in_own_pmd_1 == 1) & (m.in_own_pmd_2 == 1)).sum()
        either = ((m.in_own_pmd_1 == 1) | (m.in_own_pmd_2 == 1)).sum()
        donor_rows.append({
            'sample': sample,
            'hap_absdiff_bins_in_ref_pmd': b.absdiff[b.ref_pmd > 0.5].mean(),
            'hap_absdiff_bins_out_ref_pmd': b.absdiff[b.ref_pmd <= 0.5].mean(),
            'hap1_minus_hap2_in_ref_pmd': (m.frac_1 - m.frac_2)[m.ref_pmd_1 == 1].mean(),
            'own_pmd_hap_jaccard': both / either if either else np.nan,
        })

    per_hap = summ.merge(pd.DataFrame(rows), on=['sample', 'hap'])
    per_hap['delta_out_minus_in'] = per_hap.meth_out_ref_pmd - per_hap.meth_in_ref_pmd
    man = pd.read_csv(MANIFEST, sep='\t', usecols=['sample_id', 'superpopulation'])
    gm = pd.read_csv(GLOBAL_METH, sep='\t', usecols=['sample_id', 'hap', 'frac_meth'])
    per_hap = (per_hap.merge(man, left_on='sample', right_on='sample_id', how='left')
               .merge(gm, on=['sample_id', 'hap'], how='left')
               .drop(columns='sample_id').rename(columns={'frac_meth': 'global_meth'}))
    per_hap.to_csv(SURVEY / 'per_hap.tsv', sep='\t', index=False)

    cols = ['global_meth', 'pmd_frac_contig', 'n_pmd', 'pmd_mean_kb', 'meth_all',
            'meth_in_ref_pmd', 'meth_out_ref_pmd', 'delta_out_minus_in']
    per_donor = (per_hap.groupby(['sample', 'superpopulation'])[cols].mean().reset_index()
                 .merge(pd.DataFrame(donor_rows), on='sample', how='left'))
    per_donor.to_csv(SURVEY / 'per_donor.tsv', sep='\t', index=False)

    stats = ['meth_in_ref_pmd', 'meth_out_ref_pmd', 'delta_out_minus_in', 'pmd_frac_contig',
             'hap_absdiff_bins_in_ref_pmd', 'hap_absdiff_bins_out_ref_pmd', 'own_pmd_hap_jaccard']
    by_sp = per_donor.groupby('superpopulation')[stats].agg(['mean', 'min', 'max'])
    by_sp.columns = [f'{a}_{b}' for a, b in by_sp.columns]
    by_sp.insert(0, 'n_donors', per_donor.groupby('superpopulation').size())
    allrow = per_donor[stats].agg(['mean', 'min', 'max']).T
    allrow = pd.DataFrame({f'{s}_{k}': [allrow.loc[s, k]] for s in stats for k in allrow.columns},
                          index=['ALL'])
    allrow.insert(0, 'n_donors', len(per_donor))
    by_sp = pd.concat([by_sp, allrow])
    by_sp.index.name = 'superpopulation'
    by_sp.to_csv(SURVEY / 'by_superpop.tsv', sep='\t')

    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 30)
    show = ['n_donors'] + [f'{s}_{k}' for s in stats[:4] for k in ('mean', 'min', 'max')]
    print(by_sp[show].round(3).T.to_string())
    print('\ndonors with delta_out_minus_in <= 0:', (per_donor.delta_out_minus_in <= 0).sum())
    print('\n10 highest global-methylation donors:')
    print(per_donor.sort_values('global_meth', ascending=False).head(10)[
        ['sample', 'superpopulation', 'global_meth', 'pmd_frac_contig', 'meth_in_ref_pmd',
         'meth_out_ref_pmd', 'delta_out_minus_in']].round(3).to_string(index=False))
    print('\nhaplotype |hap1-hap2| in 10kb bins, in vs out of reference PMDs (all donors):')
    print(per_donor[stats[4:]].describe().round(4).to_string())
    sys.stdout.flush()


if __name__ == '__main__':
    main()
