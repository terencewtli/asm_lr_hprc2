#!/usr/bin/env python3
"""B01b_bin_matrix_pca.py -- genome-wide sample x 10kb-bin methylation matrix, PCA, covariate
association, per-donor continuum metrics, and a population frequency map of low-methylation
domains. Aggregates B01a's per-hap outputs. Planned analyses A/B/C in JOURNAL.md.

Bins kept for the matrix: autosomal 10kb bins with >= MIN_CPG CpGs in >= MIN_PRESENT of haps.
Donor-level value = mean of hap1/hap2 (donors missing a hap use the one they have).

A. PCA (bins centered, not scaled) on donor x bin mean methylation, two versions:
     raw        -- includes the global methylation level
     centered   -- each donor's own genome-wide bin mean subtracted first (shape only)
   Each PC is regressed on covariates one at a time (R^2): chemistry (HPRC2 Supp S6),
   superpopulation, sex, passage, where the line was established (S15), HPRC2 methylation-QC
   flag, global methylation, mean depth, ONT read N50, PMD burden.
B. Per-donor continuum metrics: global meth, fraction of bins with meth < LOW_ABS, genome-wide
   PMD burden, meth in constitutive domain bins, domain contrast (non-domain minus domain).
C. Domain frequency per bin, two definitions:
     pmd  -- frac_pmd >= 0.5 in that donor's own genome-wide dnmtools calls
     rel  -- bin meth < (donor's median bin meth - REL_DROP)   (threshold-free of the caller)
   classes: constitutive >= 0.9, variable 0.2-0.8, rare (0, 0.2), plus per-chemistry frequency.
Also a 50kb version of the PCA (10kb bins aggregated, CpG-weighted).

Outputs (results/meth_bins/): bins_index.tsv.gz, donor_bin_meth_10kb.npz, donors.tsv,
pca_{raw,centered}_{10kb,50kb}_scores.tsv, pca_*_variance.tsv, pca_covariate_r2.tsv,
donor_metrics.tsv, domain_frequency_10kb.tsv.gz (incl. fibroblast_pmd bin flag),
domain_class_summary.tsv, pmd_bin_jaccard_donors.tsv.gz (own-PMD bins, donor x donor)

Usage: B01b_bin_matrix_pca.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
BINS = PROJDIR / 'results' / 'meth_bins'
PER_HAP = BINS / 'per_hap'
MANIFEST = PROJDIR / 'tsv' / 'meta' / 'hprc2_sample_manifest.tsv'
SEQ_QC = PROJDIR / 'results' / 'qc' / 'data' / 'supp_seq_qc.csv'
COVS = PROJDIR / 'results' / 'qc' / 'data' / 'global_methylation_covariates_full_cohort.tsv'
SUPP = PROJDIR / 'reference' / 'hprc2' / 'hprc2_supp.xlsx'
FIB = PROJDIR / 'reference' / 'igvf_pgp' / 'start_merged.filt_mcg.sorted.bed.gz'
CHROMSIZES = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosomal.chromsizes')
HPRC2_METH_QC_FLAGGED = {'NA20762', 'NA19338', 'HG02583'}  # preprint Methods, Nanopore Methylation QC
MIN_CPG = 10
MIN_PRESENT = 0.95
N_PCS = 10
LOW_ABS = 0.5
REL_DROP = 0.10


def log(msg: str) -> None:
    print(msg, flush=True)


def load_matrix() -> tuple:
    files = sorted(PER_HAP.glob('*.bins10kb.tsv.gz'))
    tags = [f.name.split('.bins10kb')[0] for f in files]
    log(f'{len(files)} per-hap bin files')
    sizes = pd.read_csv(CHROMSIZES, sep='\t', header=None, names=['chrom', 'len'])
    sizes = sizes[sizes.chrom.isin([f'chr{i}' for i in range(1, 23)])]
    idx = pd.concat([pd.DataFrame({'chrom': c, 'bin_start': np.arange(0, n, 10_000)})
                     for c, n in zip(sizes.chrom, sizes.len)], ignore_index=True)
    offset = dict(zip(sizes.chrom, np.r_[0, np.cumsum((sizes.len.values + 9_999) // 10_000)[:-1]]))
    shape = (len(tags), len(idx))
    meth = np.full(shape, np.nan, dtype=np.float32)
    ncpg = np.zeros(shape, dtype=np.int32)
    pmd = np.full(shape, np.nan, dtype=np.float32)
    for i, (tag, f) in enumerate(zip(tags, files)):
        d = pd.read_csv(f, sep='\t', usecols=['chrom', 'bin_start', 'n_cpg', 'mean_meth', 'frac_pmd'])
        pos = d.chrom.map(offset).values + d.bin_start.values // 10_000
        meth[i, pos] = d.mean_meth.values
        ncpg[i, pos] = d.n_cpg.values
        pmd[i, pos] = d.frac_pmd.values
    return tags, idx, meth, ncpg, pmd


def to_donor(tags: list, x: np.ndarray) -> tuple:
    samples = [t.rsplit('_hap', 1)[0] for t in tags]
    donors = sorted(set(samples))
    out = np.full((len(donors), x.shape[1]), np.nan, dtype=np.float32)
    for j, d in enumerate(donors):
        rows = [i for i, s in enumerate(samples) if s == d]
        with np.errstate(all='ignore'):
            out[j] = np.nanmean(x[rows], axis=0)
    return donors, out


def pca(x: np.ndarray, n: int) -> tuple:
    xc = x - x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(xc, full_matrices=False)
    var = s ** 2 / (s ** 2).sum()
    n = min(n, len(s))
    return u[:, :n] * s[:n], var[:n]


def covariates(donors: list) -> pd.DataFrame:
    cov = pd.DataFrame({'sample': donors})
    man = pd.read_csv(MANIFEST, sep='\t', usecols=['sample_id', 'superpopulation'])
    seq = pd.read_csv(SEQ_QC, usecols=['sample_id', 'sequencing_chemistry_ont', 'read_N50_ont'])
    sx = pd.read_csv(COVS, sep='\t').groupby('sample_id')[['sex', 'passage']].first().reset_index()
    s15 = pd.read_excel(SUPP, sheet_name='S15', header=1)[
        ['Cell Line Sample ID', 'Location Cell Line Established']]
    s15.columns = ['sample_id', 'established']
    for t in (man, seq, sx, s15):
        cov = cov.merge(t, left_on='sample', right_on='sample_id', how='left').drop(columns='sample_id')
    cov = cov.rename(columns={'sequencing_chemistry_ont': 'chemistry'})
    cov['hprc2_meth_qc_flag'] = cov['sample'].isin(HPRC2_METH_QC_FLAGGED).astype(int)
    summ = pd.concat([pd.read_csv(f, sep='\t') for f in PER_HAP.glob('*.summary.tsv')])
    agg = summ.groupby('sample').agg(global_meth=('global_meth', 'mean'),
                                     mean_depth=('mean_depth', 'mean'),
                                     pmd_burden_bp=('pmd_native_bp', 'mean'),
                                     frac_cpg_in_pmd=('frac_cpg_in_pmd_native', 'mean'),
                                     map_rate=('map_rate', 'mean')).reset_index()
    return cov.merge(agg, on='sample', how='left')


def covariate_r2(scores: pd.DataFrame, cov: pd.DataFrame, label: str) -> pd.DataFrame:
    terms = {'chemistry': 'C(chemistry)', 'superpopulation': 'C(superpopulation)', 'sex': 'C(sex)',
             'passage': 'passage', 'established': 'C(established)',
             'hprc2_meth_qc_flag': 'hprc2_meth_qc_flag', 'global_meth': 'global_meth',
             'mean_depth': 'mean_depth', 'read_N50_ont': 'read_N50_ont',
             'frac_cpg_in_pmd': 'frac_cpg_in_pmd'}
    d = scores.merge(cov, on='sample')
    rows = []
    for pc in [c for c in scores.columns if c.startswith('PC')]:
        for name, term in terms.items():
            sub = d.dropna(subset=[name])
            if sub[name].nunique() < 2:
                continue
            r = smf.ols(f'{pc} ~ {term}', sub).fit()
            rows.append({'analysis': label, 'pc': pc, 'covariate': name, 'n': len(sub),
                         'r2': r.rsquared, 'p': r.f_pvalue})
    return pd.DataFrame(rows)


def aggregate_50kb(idx: pd.DataFrame, meth: np.ndarray, ncpg: np.ndarray) -> np.ndarray:
    grp = (idx.chrom + ':' + (idx.bin_start // 50_000).astype(str)).values
    codes, _ = pd.factorize(grp)
    w = np.where(np.isnan(meth), 0, ncpg).astype(np.float64)
    num = np.zeros((meth.shape[0], codes.max() + 1))
    den = np.zeros_like(num)
    for i in range(meth.shape[0]):
        num[i] = np.bincount(codes, weights=np.nan_to_num(meth[i]) * w[i], minlength=num.shape[1])
        den[i] = np.bincount(codes, weights=w[i], minlength=num.shape[1])
    with np.errstate(all='ignore'):
        return (num / den).astype(np.float32)


def main() -> None:
    tags, idx, meth, ncpg, pmd = load_matrix()
    ok = (ncpg >= MIN_CPG).mean(axis=0) >= MIN_PRESENT
    log(f'{len(idx):,} bins seen; {ok.sum():,} pass n_cpg>={MIN_CPG} in >={MIN_PRESENT:.0%} of haps')
    idx = idx[ok].reset_index(drop=True)
    meth, ncpg, pmd = meth[:, ok], ncpg[:, ok], pmd[:, ok]
    meth[ncpg < MIN_CPG] = np.nan

    donors, dmeth = to_donor(tags, meth)
    _, dpmd = to_donor(tags, pmd)
    _, dncpg = to_donor(tags, ncpg.astype(np.float32))
    cov = covariates(donors)
    cov.to_csv(BINS / 'donors.tsv', sep='\t', index=False)
    idx.to_csv(BINS / 'bins_index.tsv.gz', sep='\t', index=False)
    np.savez_compressed(BINS / 'donor_bin_meth_10kb.npz', meth=dmeth, pmd=dpmd,
                        donors=np.array(donors))
    log(f'donor matrix {dmeth.shape}')

    # impute the few remaining NaNs with the bin mean for PCA
    fill = np.where(np.isnan(dmeth), np.nanmean(dmeth, axis=0, keepdims=True), dmeth)
    d50 = aggregate_50kb(idx, dmeth, dncpg)
    d50 = d50[:, ~np.isnan(d50).any(axis=0)]
    r2_all = []
    for res, x in (('10kb', fill), ('50kb', d50)):
        for label, xx in (('raw', x), ('centered', x - x.mean(axis=1, keepdims=True))):
            sc, var = pca(xx, N_PCS)
            scores = pd.DataFrame(sc, columns=[f'PC{i + 1}' for i in range(sc.shape[1])])
            scores.insert(0, 'sample', donors)
            tag = f'{label}_{res}'
            scores.to_csv(BINS / f'pca_{tag}_scores.tsv', sep='\t', index=False)
            pd.DataFrame({'pc': scores.columns[1:], 'var_explained': var}).to_csv(
                BINS / f'pca_{tag}_variance.tsv', sep='\t', index=False)
            log(f'PCA {tag}: variance explained ' + ' '.join(f'{v:.3f}' for v in var))
            r2_all.append(covariate_r2(scores, cov, tag))
    r2 = pd.concat(r2_all)
    r2.to_csv(BINS / 'pca_covariate_r2.tsv', sep='\t', index=False)
    for tag in ('raw_10kb', 'centered_10kb'):
        wide = r2[(r2.analysis == tag) & r2.pc.isin([f'PC{i}' for i in range(1, 6)])].pivot(
            index='covariate', columns='pc', values='r2')
        log(f'\nPC ~ covariate R^2 ({tag}):\n' + wide.round(3).to_string())

    # C. domain frequency
    donor_median = np.nanmedian(dmeth, axis=1, keepdims=True)
    low_pmd = np.where(np.isnan(dpmd), np.nan, (dpmd >= 0.5).astype(np.float32))
    low_rel = np.where(np.isnan(dmeth), np.nan, (dmeth < donor_median - REL_DROP).astype(np.float32))
    freq = idx.copy()
    freq['n_donors'] = (~np.isnan(dmeth)).sum(axis=0)
    freq['mean_meth'] = np.nanmean(dmeth, axis=0)
    freq['sd_meth'] = np.nanstd(dmeth, axis=0)
    freq['freq_pmd'] = np.nanmean(low_pmd, axis=0)
    freq['freq_rel'] = np.nanmean(low_rel, axis=0)
    chem = cov.chemistry.values
    for c in ('R941', 'R1041'):
        sel = chem == c
        freq[f'freq_pmd_{c}'] = np.nanmean(low_pmd[sel], axis=0)
        freq[f'freq_rel_{c}'] = np.nanmean(low_rel[sel], axis=0)

    def classify(f: pd.Series) -> pd.Series:
        return pd.cut(f, [-0.001, 0, 0.2, 0.8, 0.9, 1.0],
                      labels=['never', 'rare', 'variable', 'common', 'constitutive'])

    freq['class_pmd'] = classify(freq.freq_pmd)
    freq['class_rel'] = classify(freq.freq_rel)
    freq.round(4).to_csv(BINS / 'domain_frequency_10kb.tsv.gz', sep='\t', index=False)
    cls = pd.concat([freq.class_pmd.value_counts(normalize=True).rename('pmd'),
                     freq.class_rel.value_counts(normalize=True).rename('rel')], axis=1)
    cls.to_csv(BINS / 'domain_class_summary.tsv', sep='\t')
    log('\ndomain class fraction of bins:\n' + cls.round(3).to_string())
    log(f'R941 vs R1041 domain frequency correlation: pmd r='
        f'{np.corrcoef(freq.freq_pmd_R941, freq.freq_pmd_R1041)[0, 1]:.3f}, rel r='
        f'{np.corrcoef(freq.freq_rel_R941, freq.freq_rel_R1041)[0, 1]:.3f}')

    # bin-level Jaccard of own-PMD bins: donor vs donor, donor vs fibroblast (igvf_pgp)
    fibb = pd.read_csv(FIB, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 's', 'e'])
    fib_cov = np.zeros(len(idx))
    key = {k: i for i, k in enumerate(zip(idx.chrom, idx.bin_start))}
    for c, s, e in zip(fibb.chrom, fibb.s, fibb.e):
        for b0 in range((s // 10_000) * 10_000, e, 10_000):
            i = key.get((c, b0))
            if i is not None:
                fib_cov[i] += min(e, b0 + 10_000) - max(s, b0)
    fib_bin = fib_cov >= 5_000
    lp = np.nan_to_num(low_pmd).astype(bool)
    inter = lp.astype(np.float32) @ lp.T.astype(np.float32)
    tot = lp.sum(axis=1)
    jac = inter / (tot[:, None] + tot[None, :] - inter)
    pd.DataFrame(jac, index=donors, columns=donors).round(4).to_csv(BINS / 'pmd_bin_jaccard_donors.tsv.gz', sep='\t')
    jf = (lp & fib_bin).sum(axis=1) / (lp | fib_bin).sum(axis=1)
    iu = np.triu_indices(len(donors), 1)
    log(f'\nPMD-bin Jaccard donor-vs-donor: median {np.median(jac[iu]):.3f} '
        f'(IQR {np.percentile(jac[iu], 25):.3f}-{np.percentile(jac[iu], 75):.3f}); '
        f'donor-vs-fibroblast: median {np.median(jf):.3f} (range {jf.min():.3f}-{jf.max():.3f}); '
        f'fibroblast PMD bins {fib_bin.mean():.3f} of kept bins')
    freq['fibroblast_pmd'] = fib_bin.astype(int)
    freq.round(4).to_csv(BINS / 'domain_frequency_10kb.tsv.gz', sep='\t', index=False)
    log('\nfibroblast PMD bin fraction by LCL domain class:\n'
        + freq.groupby('class_pmd').fibroblast_pmd.mean().round(3).to_string())

    # B. per-donor continuum metrics
    const = (freq.class_pmd == 'constitutive').values
    never = (freq.class_pmd == 'never').values
    met = cov[['sample', 'chemistry', 'superpopulation', 'global_meth', 'frac_cpg_in_pmd',
               'pmd_burden_bp']].copy()
    met['frac_bins_below_0.5'] = np.nanmean(dmeth < LOW_ABS, axis=1)
    met['frac_bins_rel_low'] = np.nanmean(low_rel, axis=1)
    met['meth_constitutive_domains'] = np.nanmean(dmeth[:, const], axis=1)
    met['meth_never_domain_bins'] = np.nanmean(dmeth[:, never], axis=1)
    met['domain_contrast'] = met.meth_never_domain_bins - met.meth_constitutive_domains
    for d in ('NA19338', 'HG04187'):
        if d in donors:
            j = donors.index(d)
            own = low_pmd[j] == 1
            log(f'{d}: own PMD bins by population class: '
                + freq.class_pmd[own].value_counts(normalize=True).round(3).to_dict().__repr__())
    met.to_csv(BINS / 'donor_metrics.tsv', sep='\t', index=False)
    num = met.select_dtypes('number')
    log('\ncontinuum metric correlations:\n' + num.corr().round(3).to_string())
    for y in ('global_meth', 'frac_cpg_in_pmd', 'domain_contrast'):
        sub = met.dropna(subset=[y, 'chemistry', 'superpopulation'])
        for f in ('C(chemistry)', 'C(superpopulation)', 'C(chemistry) + C(superpopulation)'):
            r = smf.ols(f'{y} ~ {f}', sub).fit()
            log(f'{y} ~ {f}: R2={r.rsquared:.3f} p={r.f_pvalue:.2g} (n={len(sub)})')


if __name__ == '__main__':
    main()
