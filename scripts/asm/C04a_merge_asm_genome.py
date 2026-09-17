#!/usr/bin/env python3
"""C04a_merge_asm_genome.py -- genome-wide merge of C02a's per-(sample, chrom) ASM calls.

Streams one donor at a time (all donors x ~2M regions is too big to hold at once):
  - BH-FDR across the donor's whole genome, separately per region set (cpg, zink)
  - ASM call: q < Q_MAX and |delta| >= MIN_DELTA
  - writes that donor's calls-only table, and adds the donor to per-region counters split by
    ONT chemistry (HPRC2 Supp S6)

Outputs (results/asm/genome/):
  per_donor/<sample>.asm_sig.tsv.gz    significant regions for that donor (all columns + q)
  donor_summary.tsv                    per donor: n_tested, n_asm, frac, median reads, chemistry,
                                       superpopulation, informative-read yield (C02a logs)
  region_penetrance.tsv.gz             per region: n_tested / n_asm overall and per chemistry,
                                       mean |delta| among ASM calls
  zink_control.tsv                     per donor x Zink DMR: delta, p, q, called
  zink_control_summary.tsv             Zink detection rate / |delta| by chemistry (+ tests)

Usage: C04a_merge_asm_genome.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
CALLS = PROJDIR / 'results' / 'asm' / 'calls'
REGIONS = PROJDIR / 'results' / 'asm' / 'regions' / 'regions_hg38.tsv.gz'
OUT = PROJDIR / 'results' / 'asm' / 'genome'
SEQ_QC = PROJDIR / 'results' / 'qc' / 'data' / 'supp_seq_qc.csv'
MANIFEST = PROJDIR / 'tsv' / 'meta' / 'hprc2_sample_manifest.tsv'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
Q_MAX = 0.05
MIN_DELTA = 0.2


def main() -> None:
    (OUT / 'per_donor').mkdir(parents=True, exist_ok=True)
    reg = pd.read_csv(REGIONS, sep='\t', usecols=['region_id', 'set', 'chrom', 'start', 'end', 'n_cpg_ref'])
    ridx = pd.Series(np.arange(len(reg)), index=reg.region_id)
    chem = pd.read_csv(SEQ_QC).set_index('sample_id').sequencing_chemistry_ont
    sp = pd.read_csv(MANIFEST, sep='\t').set_index('sample_id').superpopulation
    groups = ['all', 'R941', 'R1041']
    n_test = {g: np.zeros(len(reg), np.int32) for g in groups}
    n_asm = {g: np.zeros(len(reg), np.int32) for g in groups}
    sum_absd = np.zeros(len(reg))

    samples = sorted({f.name.split('_chr')[0] for f in CALLS.glob('*.asm.tsv.gz')})
    donor_rows, zink_rows = [], []
    for s in samples:
        files = [CALLS / f'{s}_{c}.asm.tsv.gz' for c in AUTOSOMES]
        missing = [f.name for f in files if not f.exists()]
        d = pd.concat([pd.read_csv(f, sep='\t') for f in files if f.exists()], ignore_index=True)
        logs = pd.concat([pd.read_csv(f.with_name(f.name.replace('.asm.tsv.gz', '.log.tsv')), sep='\t')
                          for f in files if f.with_name(f.name.replace('.asm.tsv.gz', '.log.tsv')).exists()])
        d['q'] = np.nan
        for rset, g in d.groupby('set'):
            d.loc[g.index, 'q'] = multipletests(g.pval.values, method='fdr_bh')[1]
        d['asm'] = (d.q < Q_MAX) & (d.delta.abs() >= MIN_DELTA)
        d[d.asm].to_csv(OUT / 'per_donor' / f'{s}.asm_sig.tsv.gz', sep='\t', index=False, float_format='%.6g')
        c = chem.get(s, 'NA')
        ii = ridx[d.region_id].values
        for g in ('all', c):
            if g in n_test:
                np.add.at(n_test[g], ii, 1)
                np.add.at(n_asm[g], ii[d.asm.values], 1)
        np.add.at(sum_absd, ii[d.asm.values], d.delta.abs().values[d.asm.values])
        cg = d[d.set == 'cpg']
        donor_rows.append({'sample': s, 'chemistry': c, 'superpopulation': sp.get(s, 'NA'),
                           'n_chrom': 22 - len(missing), 'n_tested': len(cg), 'n_asm': int(cg.asm.sum()),
                           'frac_asm': cg.asm.mean(), 'median_reads': float(np.median(cg.n_reads1 + cg.n_reads2)),
                           'frac_reads_het': logs.reads_het.sum() / max(logs.reads.sum(), 1),
                           'frac_reads_kept': logs.reads_kept.sum() / max(logs.reads.sum(), 1),
                           'frac_calls_cpg': logs.calls_cpg.sum() / max(logs.calls.sum(), 1)})
        z = d[d.set == 'zink'][['region_id', 'n_reads1', 'n_reads2', 'delta', 'pval', 'q', 'asm']]
        zink_rows.append(z.assign(sample=s, chemistry=c))
        print(f'{s}: {len(cg):,} tested, {int(cg.asm.sum()):,} ASM; missing {len(missing)} chroms', flush=True)

    ds = pd.DataFrame(donor_rows)
    ds.to_csv(OUT / 'donor_summary.tsv', sep='\t', index=False)
    rp = reg.copy()
    for g in groups:
        rp[f'n_tested_{g}'] = n_test[g]
        rp[f'n_asm_{g}'] = n_asm[g]
        with np.errstate(all='ignore'):
            rp[f'frac_asm_{g}'] = n_asm[g] / n_test[g]
    with np.errstate(all='ignore'):
        rp['mean_absdelta_asm'] = sum_absd / n_asm['all']
    rp[rp.n_tested_all > 0].to_csv(OUT / 'region_penetrance.tsv.gz', sep='\t', index=False, float_format='%.4g')

    z = pd.concat(zink_rows, ignore_index=True)
    z.to_csv(OUT / 'zink_control.tsv', sep='\t', index=False)
    zd = z.groupby(['sample', 'chemistry']).agg(n=('asm', 'size'), det=('asm', 'mean'),
                                                absd=('delta', lambda x: x.abs().mean())).reset_index()
    rows = []
    for metric in ('det', 'absd'):
        a, b = zd[zd.chemistry == 'R941'][metric], zd[zd.chemistry == 'R1041'][metric]
        p = mannwhitneyu(a, b).pvalue if len(a) and len(b) else np.nan
        rows.append({'metric': metric, 'R941_mean': a.mean(), 'R941_n': len(a),
                     'R1041_mean': b.mean(), 'R1041_n': len(b), 'mannwhitney_p': p})
    zs = pd.DataFrame(rows)
    zs.to_csv(OUT / 'zink_control_summary.tsv', sep='\t', index=False)
    print('\nZink imprinted-DMR control by chemistry (det = fraction of DMRs called ASM per donor, '
          'absd = mean |hap1-hap2|):\n' + zs.round(4).to_string(index=False))
    print('\ndonor ASM yield by chemistry:\n' + ds.groupby('chemistry')[
        ['n_tested', 'n_asm', 'frac_asm', 'median_reads', 'frac_reads_het']].median().round(4).to_string())
    print('\ndonor ASM yield by superpopulation:\n' + ds.groupby('superpopulation')[
        ['n_tested', 'n_asm', 'frac_asm', 'frac_reads_het']].median().round(4).to_string())


if __name__ == '__main__':
    main()
