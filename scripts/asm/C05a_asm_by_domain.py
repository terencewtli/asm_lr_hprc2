#!/usr/bin/env python3
"""C05a_asm_by_domain.py -- is allele-specific methylation concentrated in PMDs, and does it look
genetic or stochastic? Runs on C04a's genome-wide merge.

Joins each tested region (C04a region_penetrance) to:
  - LCL domain class and mean mCG (B01b domain_frequency_10kb)
  - replication timing / LAD (B06a rt_lad_10kb)
  - HPRC2 promoter mQTL bins (reference/hprc2/s11_promoter_mqtl_hg38.bed)
and reports:
  1. ASM rate (fraction of tested donors called ASM) and |delta| by domain class and RT decile
  2. recurrence: regions ASM in many donors (consistent with a cis-genetic cause, since donors
     share variants) vs donor-private regions (consistent with stochastic/epigenetic drift),
     per domain class
  3. whether recurrent-ASM regions are enriched for HPRC2 promoter mQTLs relative to tested
     regions matched on domain class and CpG count
  4. per donor: ASM count vs domain depth and chemistry (does a deeper domain state produce more
     ASM calls, i.e. is ASM yield confounded by the methylome state?)

Outputs (results/asm/genome/): asm_by_domain_class.tsv, asm_by_rt_decile.tsv,
asm_recurrence_by_class.tsv, asm_mqtl_enrichment.tsv, asm_per_donor_vs_depth.tsv

Usage: C05a_asm_by_domain.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
GEN = PROJDIR / 'results' / 'asm' / 'genome'
BINS = PROJDIR / 'results' / 'meth_bins'
MQTL = PROJDIR / 'reference' / 'hprc2' / 's11_promoter_mqtl_hg38.bed'
QC13 = PROJDIR / 'results' / 'qc' / 'data' / 'qc13' / 'donor_pmd_expansion_metrics.tsv'
BIN = 10_000
CLASSES = ['never', 'rare', 'variable', 'common', 'constitutive']


def log(msg: str) -> None:
    print(msg, flush=True)


def main() -> None:
    rp = pd.read_csv(GEN / 'region_penetrance.tsv.gz', sep='\t')
    rp = rp[rp['set'] == 'cpg'].copy()
    rp['bin_start'] = (rp.start // BIN) * BIN
    freq = pd.read_csv(BINS / 'domain_frequency_10kb.tsv.gz', sep='\t',
                       usecols=['chrom', 'bin_start', 'class_pmd', 'mean_meth', 'freq_pmd'])
    rt = pd.read_csv(BINS / 'annotations' / 'rt_lad_10kb.tsv.gz', sep='\t')
    mq = pd.read_csv(MQTL, sep='\t', header=None, usecols=[0, 1], names=['chrom', 'pos'])
    mq['bin_start'] = (mq.pos // BIN) * BIN
    mq_bins = set(map(tuple, mq[['chrom', 'bin_start']].drop_duplicates().values))

    d = rp.merge(freq, on=['chrom', 'bin_start'], how='left').merge(rt, on=['chrom', 'bin_start'], how='left')
    d['has_mqtl'] = [(c, b) in mq_bins for c, b in zip(d.chrom, d.bin_start)]
    d['asm_rate'] = d.n_asm_all / d.n_tested_all
    d = d[d.n_tested_all >= 50]
    log(f'{len(d):,} regions tested in >=50 donors; {d.n_asm_all.sum():,} donor-level ASM calls')

    by_class = d.groupby('class_pmd').agg(
        n_regions=('asm_rate', 'size'), median_tested=('n_tested_all', 'median'),
        mean_asm_rate=('asm_rate', 'mean'), frac_regions_any_asm=('n_asm_all', lambda x: (x > 0).mean()),
        mean_absdelta=('mean_absdelta_asm', 'mean'), n_cpg=('n_cpg_ref', 'median')).reindex(CLASSES)
    by_class.to_csv(GEN / 'asm_by_domain_class.tsv', sep='\t')
    log('\nASM by domain class:\n' + by_class.round(4).to_string())

    d['rt_decile'] = pd.qcut(d.rt, 10, labels=False, duplicates='drop')
    by_rt = d.groupby('rt_decile').agg(n=('asm_rate', 'size'), rt=('rt', 'median'),
                                       mean_asm_rate=('asm_rate', 'mean'),
                                       frac_constitutive=('class_pmd', lambda x: (x == 'constitutive').mean()))
    by_rt.to_csv(GEN / 'asm_by_rt_decile.tsv', sep='\t')
    log('\nASM by replication-timing decile (low = late):\n' + by_rt.round(4).to_string())

    # recurrence
    d['recurrence'] = pd.cut(d.asm_rate, [-0.001, 0, 0.02, 0.1, 0.5, 1.0],
                             labels=['none', 'private (<2%)', 'uncommon (2-10%)', 'common (10-50%)', 'shared (>50%)'])
    rec = pd.crosstab(d.class_pmd, d.recurrence, normalize='index').reindex(CLASSES)
    rec.to_csv(GEN / 'asm_recurrence_by_class.tsv', sep='\t')
    log('\nfraction of regions by ASM recurrence, per domain class:\n' + rec.round(4).to_string())

    # mQTL enrichment of recurrent ASM, matched on class and CpG count
    d['cpg_bin'] = pd.qcut(d.n_cpg_ref, 5, labels=False, duplicates='drop')
    d['recurrent'] = d.asm_rate >= 0.1
    fit = smf.logit('has_mqtl ~ recurrent + C(class_pmd) + C(cpg_bin)', d.dropna(subset=['class_pmd'])).fit(disp=0)
    enr = pd.DataFrame({'term': fit.params.index, 'OR': np.exp(fit.params.values), 'p': fit.pvalues.values})
    enr.to_csv(GEN / 'asm_mqtl_enrichment.tsv', sep='\t', index=False)
    log('\nlogit: has HPRC2 promoter mQTL ~ recurrent ASM + class + CpG quintile\n'
        + enr[enr.term.str.contains('recurrent|Intercept')].round(4).to_string(index=False))
    log('\nraw mQTL-bin fraction by recurrence:\n'
        + d.groupby('recurrence').has_mqtl.agg(['size', 'mean']).round(4).to_string())

    # per donor
    ds = pd.read_csv(GEN / 'donor_summary.tsv', sep='\t')
    dm = pd.read_csv(QC13, sep='\t')[['sample', 'depth_constitutive', 'breadth_rel', 'mcg_never']]
    pd_ = ds.merge(dm, on='sample', how='left')
    pd_.to_csv(GEN / 'asm_per_donor_vs_depth.tsv', sep='\t', index=False)
    sub = pd_.dropna(subset=['depth_constitutive', 'frac_asm'])
    if len(sub) > 10:
        f1 = smf.ols('frac_asm ~ depth_constitutive + C(chemistry) + median_reads', sub).fit()
        log('\nper-donor ASM fraction ~ domain depth + chemistry + depth of coverage:\n'
            + pd.DataFrame({'coef': f1.params, 'p': f1.pvalues}).round(5).to_string())
        log(f'R2 {f1.rsquared:.3f}; corr(frac_asm, depth_constitutive) '
            f'{np.corrcoef(sub.frac_asm, sub.depth_constitutive)[0, 1]:.3f}')


if __name__ == '__main__':
    main()
