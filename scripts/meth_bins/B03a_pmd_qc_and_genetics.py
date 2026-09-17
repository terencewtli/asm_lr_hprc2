#!/usr/bin/env python3
"""B03a_pmd_qc_and_genetics.py -- genome-wide PMD spot-check + planned analysis D (genetically
regulated methylation vs PMD/domain class). Uses B01a/B01b outputs only.

1. PMD QC per haplotype (B01a summaries + native PMD beds): burden, size distribution, mCG inside
   vs outside, by chemistry; hap1 vs hap2 agreement per donor (Jaccard of own-PMD 10kb bins;
   mean |hap1-hap2| of 10kb mCG in constitutive vs never-PMD bins).
2. D(a) HPRC2 promoter mQTLs (Supp S11, lifted CHM13->hg38) by LCL domain class: among 10kb bins
   containing a protein-coding/lncRNA TSS, the fraction with >=1 mQTL bin, per class; logistic
   regression of has_mQTL on class adjusted for log CpG count and n TSS.
3. D(b) per-bin variance decomposition across donors (10kb mean mCG, haps kept separate):
     between-donor variance of donor means, within-donor variance (hap1-hap2)^2/2,
     R^2 of donor means on donor global mCG ("state") and on state + chemistry;
   summarized by domain class. Predicted: domain bins dominated by state; never-PMD bins with
   relatively more within-donor (cis/allelic) variance.

Outputs (results/meth_bins/qc_genetics/): pmd_qc_per_hap.tsv, pmd_size_quantiles.tsv,
hap_agreement_per_donor.tsv, mqtl_by_domain_class.tsv, mqtl_logit.txt,
variance_by_domain_class.tsv, bin_variance_10kb.tsv.gz

Usage: B03a_pmd_qc_and_genetics.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
BINS = PROJDIR / 'results' / 'meth_bins'
PER_HAP = BINS / 'per_hap'
PMD_DIR = PROJDIR / 'data' / 'pmds'
OUT = BINS / 'qc_genetics'
SEQ_QC = PROJDIR / 'results' / 'qc' / 'data' / 'supp_seq_qc.csv'
MQTL = PROJDIR / 'reference' / 'hprc2' / 's11_promoter_mqtl_hg38.bed'
TSS = PROJDIR / 'reference' / 'hprc2' / 'gencode_v43_autosome_gene_tss.bed'
BIN = 10_000
CLASSES = ['never', 'rare', 'variable', 'common', 'constitutive']


def log(msg: str) -> None:
    print(msg, flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    chem = pd.read_csv(SEQ_QC).set_index('sample_id').sequencing_chemistry_ont
    freq = pd.read_csv(BINS / 'domain_frequency_10kb.tsv.gz', sep='\t')
    key = pd.MultiIndex.from_frame(freq[['chrom', 'bin_start']])

    # 1. PMD QC
    summ = pd.concat([pd.read_csv(f, sep='\t') for f in PER_HAP.glob('*.summary.tsv')])
    summ['chemistry'] = summ['sample'].map(chem)
    summ['contrast_native'] = summ.meth_out_pmd_native - summ.meth_in_pmd_native
    summ.to_csv(OUT / 'pmd_qc_per_hap.tsv', sep='\t', index=False)
    cols = ['pmd_native_bp', 'frac_cpg_in_pmd_native', 'meth_in_pmd_native', 'meth_out_pmd_native',
            'contrast_native', 'global_wmeth', 'mean_depth', 'map_rate']
    log('PMD QC per haplotype, median [min, max]:')
    for c in cols:
        log(f'  {c}: {summ[c].median():.4g} [{summ[c].min():.4g}, {summ[c].max():.4g}]  '
            f'R941 {summ[summ.chemistry == "R941"][c].median():.4g}  R1041 {summ[summ.chemistry == "R1041"][c].median():.4g}')
    sizes = []
    for f in PMD_DIR.glob('*/*.pmd.bed'):
        if f.stat().st_size:
            b = pd.read_csv(f, sep='\t', header=None, usecols=[1, 2])
            sizes.append((b[2] - b[1]).values)
    sizes = np.concatenate(sizes)
    q = pd.Series(np.quantile(sizes, [.05, .25, .5, .75, .95, .99]), index=['q05', 'q25', 'q50', 'q75', 'q95', 'q99'])
    q['n_domains'] = len(sizes)
    q['frac_bp_in_domains_gt_1Mb'] = sizes[sizes > 1e6].sum() / sizes.sum()
    q.to_csv(OUT / 'pmd_size_quantiles.tsv', sep='\t', header=['value'])
    log('PMD size (bp) quantiles pooled over haplotypes:\n' + q.round(3).to_string())

    # per-hap 10kb matrices on the B01b bin index
    files = sorted(PER_HAP.glob('*.bins10kb.tsv.gz'))
    tags = [f.name[:-len('.bins10kb.tsv.gz')] for f in files]
    M = np.full((len(tags), len(freq)), np.nan, dtype=np.float32)
    PM = np.full_like(M, np.nan)
    NC_sum = np.zeros(len(freq))
    NC_n = np.zeros(len(freq))
    for i, f in enumerate(files):
        d = pd.read_csv(f, sep='\t', usecols=['chrom', 'bin_start', 'n_cpg', 'mean_meth', 'frac_pmd'])
        d = d[d.n_cpg >= 10]
        j = key.get_indexer(pd.MultiIndex.from_frame(d[['chrom', 'bin_start']]))
        ok = j >= 0
        M[i, j[ok]] = d.mean_meth.values[ok]
        PM[i, j[ok]] = d.frac_pmd.values[ok]
        NC_sum[j[ok]] += d.n_cpg.values[ok]
        NC_n[j[ok]] += 1
    samples = np.array([t.rsplit('_hap', 1)[0] for t in tags])
    haps = np.array([int(t.rsplit('_hap', 1)[1]) for t in tags])
    donors = sorted(set(samples[haps == 1]) & set(samples[haps == 2]))
    i1 = np.array([tags.index(f'{d}_hap1') for d in donors])
    i2 = np.array([tags.index(f'{d}_hap2') for d in donors])
    cls = freq.class_pmd.values
    rows = []
    for d, a, b in zip(donors, i1, i2):
        pa, pb = PM[a] >= 0.5, PM[b] >= 0.5
        valid = ~np.isnan(PM[a]) & ~np.isnan(PM[b])
        diff = np.abs(M[a] - M[b])
        rows.append({'sample': d, 'chemistry': chem.get(d, 'NA'),
                     'pmd_bin_jaccard_hap1_hap2': (pa & pb & valid).sum() / max((((pa | pb) & valid)).sum(), 1),
                     'absdiff_constitutive': np.nanmean(diff[cls == 'constitutive']),
                     'absdiff_never': np.nanmean(diff[cls == 'never']),
                     'absdiff_variable': np.nanmean(diff[cls == 'variable'])})
    ha = pd.DataFrame(rows)
    ha.to_csv(OUT / 'hap_agreement_per_donor.tsv', sep='\t', index=False)
    log('\nhap1 vs hap2 per donor (median [min, max]):\n' + ha.describe().loc[['50%', 'min', 'max']].round(4).to_string())

    # 2. mQTL enrichment by domain class
    mq = pd.read_csv(MQTL, sep='\t', header=None, usecols=[0, 1], names=['chrom', 'start'])
    mq['bin_start'] = (mq.start // BIN) * BIN
    mqn = mq.groupby(['chrom', 'bin_start']).size()
    tss = pd.read_csv(TSS, sep='\t', header=None, usecols=[0, 1, 4], names=['chrom', 'pos', 'type'])
    tss = tss[tss.type.isin(['protein_coding', 'lncRNA'])]
    tss['bin_start'] = (tss.pos // BIN) * BIN
    tn = tss.groupby(['chrom', 'bin_start']).size()
    fb = freq[['chrom', 'bin_start', 'class_pmd', 'class_rel', 'mean_meth', 'sd_meth']].copy()
    fb['n_mqtl'] = mqn.reindex(key).fillna(0).values
    fb['n_tss'] = tn.reindex(key).fillna(0).values
    n_cpg_bin = NC_sum / np.maximum(NC_n, 1)
    fb['log_ncpg'] = np.log10(n_cpg_bin + 1)
    fb['has_mqtl'] = (fb.n_mqtl > 0).astype(int)
    t = fb[fb.n_tss > 0].copy()
    tab = t.groupby('class_pmd').agg(n_bins=('has_mqtl', 'size'), frac_with_mqtl=('has_mqtl', 'mean'),
                                     mqtl_per_bin=('n_mqtl', 'mean'), median_log_ncpg=('log_ncpg', 'median')).reindex(CLASSES)
    tab.to_csv(OUT / 'mqtl_by_domain_class.tsv', sep='\t')
    log(f'\nHPRC2 promoter mQTL bins by LCL domain class (TSS-containing 10kb bins, n={len(t):,}; '
        f'{int(fb.n_mqtl.sum()):,} mQTL bins placed):\n' + tab.round(4).to_string())
    t['class_pmd'] = pd.Categorical(t.class_pmd, categories=CLASSES)
    fit = smf.logit('has_mqtl ~ C(class_pmd) + log_ncpg + np.log1p(n_tss)', t).fit(disp=0)
    (OUT / 'mqtl_logit.txt').write_text(fit.summary().as_text())
    log('\nlogit has_mqtl (ref = never), odds ratios:\n'
        + pd.DataFrame({'OR': np.exp(fit.params), 'p': fit.pvalues}).round(4).to_string())

    # 3. variance decomposition
    D1, D2 = M[i1], M[i2]
    dm = (D1 + D2) / 2
    v_between = np.nanvar(dm, axis=0, ddof=1)
    v_within = np.nanmean((D1 - D2) ** 2, axis=0) / 2
    gm = np.nanmean(dm, axis=1)
    ch = np.array([1.0 if chem.get(d) == 'R1041' else 0.0 for d in donors])
    X1 = np.column_stack([np.ones(len(donors)), gm])
    X2 = np.column_stack([X1, ch])
    ok_bins = np.isnan(dm).sum(axis=0) == 0
    Y = dm[:, ok_bins]
    Yc = Y - Y.mean(axis=0)
    sst = (Yc ** 2).sum(axis=0)
    r2 = {}
    for name, X in (('state', X1), ('state_chem', X2), ('chem', np.column_stack([np.ones(len(donors)), ch]))):
        beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
        r2[name] = 1 - ((Y - X @ beta) ** 2).sum(axis=0) / sst
    bv = fb[['chrom', 'bin_start', 'class_pmd', 'class_rel', 'n_mqtl', 'n_tss']].copy()
    bv['v_between'] = v_between
    bv['v_within'] = v_within
    for k, v in r2.items():
        bv[f'r2_{k}'] = np.nan
        bv.loc[ok_bins, f'r2_{k}'] = v
    bv.round(5).to_csv(OUT / 'bin_variance_10kb.tsv.gz', sep='\t', index=False)
    vt = bv.groupby('class_pmd').agg(n=('v_between', 'size'), v_between=('v_between', 'median'),
                                     v_within=('v_within', 'median'), r2_state=('r2_state', 'median'),
                                     r2_chem=('r2_chem', 'median'), r2_state_chem=('r2_state_chem', 'median')).reindex(CLASSES)
    vt['within_over_between'] = vt.v_within / vt.v_between
    vt.to_csv(OUT / 'variance_by_domain_class.tsv', sep='\t')
    log(f'\nper-bin variance across {len(donors)} donors, medians by domain class:\n' + vt.round(5).to_string())
    mt = bv[bv.n_tss > 0].assign(mqtl=lambda x: x.n_mqtl > 0).groupby(['class_pmd', 'mqtl']).agg(
        n=('v_between', 'size'), v_between=('v_between', 'median'), v_within=('v_within', 'median'),
        r2_state=('r2_state', 'median')).reindex(CLASSES, level=0)
    log('\nTSS bins, with vs without an HPRC2 mQTL:\n' + mt.round(5).to_string())


if __name__ == '__main__':
    main()
