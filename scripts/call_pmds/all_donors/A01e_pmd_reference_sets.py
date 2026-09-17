#!/usr/bin/env python3
"""A01e_pmd_reference_sets.py -- hg38 PMD reference interval sets shared by every donor, for
metagene profiles (A01f) and LCL-vs-fibroblast comparisons.

Sets (results/pmd_metagene/ref_sets/<name>.bed, chrom/start/end, sorted):
  lcl_constitutive     10kb bins with class_pmd == constitutive (own genome-wide PMD in >=90% of
                       donors, B01b), adjacent bins merged allowing one missing/other bin gap
  lcl_common_plus      same with freq_pmd >= 0.8
  fibroblast           igvf_pgp merged fibroblast (start) PMDs
  NA19338_consensus    NA19338 hap1 AND hap2 hg38 PMDs (B01a), merged within 5kb
  HG04187_consensus    same for the highest-global-mCG donor
  NA20762_consensus    same for the lowest-global-mCG donor
Each set is also written as <name>_ge300kb.bed (intervals >= MIN_SIZE).

Usage: A01e_pmd_reference_sets.py
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
BINS = PROJDIR / 'results' / 'meth_bins'
PER_HAP = BINS / 'per_hap'
OUT = PROJDIR / 'results' / 'pmd_metagene' / 'ref_sets'
FIB = PROJDIR / 'reference' / 'igvf_pgp' / 'start_merged.filt_mcg.sorted.bed.gz'
BEDTOOLS = '/u/home/t/terencew/bin/bedtools'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
BIN = 10_000
MIN_SIZE = 300_000
DONORS = ['NA19338', 'HG04187', 'NA20762']


def bins_to_intervals(fr: pd.DataFrame, mask: np.ndarray) -> pd.DataFrame:
    out = []
    for c, g in fr[mask].groupby('chrom', sort=False):
        s = g.bin_start.values
        brk = np.r_[True, np.diff(s) > 2 * BIN]       # tolerate one missing bin
        gid = np.cumsum(brk)
        d = pd.DataFrame({'s': s, 'g': gid}).groupby('g').s.agg(['min', 'max'])
        out.append(pd.DataFrame({'chrom': c, 'start': d['min'].values, 'end': d['max'].values + BIN}))
    return pd.concat(out, ignore_index=True)


def write(df: pd.DataFrame, name: str) -> None:
    df = df[df.chrom.isin(AUTOSOMES)].copy()
    df['o'] = df.chrom.str[3:].astype(int)
    df = df.sort_values(['o', 'start']).drop(columns='o')
    df[['chrom', 'start', 'end']].to_csv(OUT / f'{name}.bed', sep='\t', header=False, index=False)
    big = df[df.end - df.start >= MIN_SIZE]
    big[['chrom', 'start', 'end']].to_csv(OUT / f'{name}_ge300kb.bed', sep='\t', header=False, index=False)
    L = df.end - df.start
    print(f'{name}: {len(df):,} intervals, {L.sum() / 1e6:.0f} Mb, median {L.median() / 1e3:.0f} kb; '
          f'>=300kb: {len(big):,} ({(big.end - big.start).sum() / 1e6:.0f} Mb)', flush=True)


def consensus(sample: str) -> pd.DataFrame:
    a, b = (PER_HAP / f'{sample}_hap{h}.hg38.pmd.bed' for h in (1, 2))
    tmp = OUT / f'.{sample}.tmp.bed'
    inter = subprocess.run([BEDTOOLS, 'intersect', '-a', str(a), '-b', str(b)], check=True,
                           capture_output=True, text=True).stdout
    tmp.write_text(inter)
    srt = subprocess.run([BEDTOOLS, 'sort', '-i', str(tmp)], check=True, capture_output=True, text=True).stdout
    tmp.write_text(srt)
    mg = subprocess.run([BEDTOOLS, 'merge', '-d', '5000', '-i', str(tmp)], check=True,
                        capture_output=True, text=True).stdout
    tmp.unlink()
    rows = [l.split('\t')[:3] for l in mg.strip().split('\n') if l]
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end']).astype({'start': int, 'end': int})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fr = pd.read_csv(BINS / 'domain_frequency_10kb.tsv.gz', sep='\t')
    write(bins_to_intervals(fr, (fr.class_pmd == 'constitutive').values), 'lcl_constitutive')
    write(bins_to_intervals(fr, (fr.freq_pmd >= 0.8).values), 'lcl_common_plus')
    fib = pd.read_csv(FIB, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    write(fib, 'fibroblast')
    for d in DONORS:
        write(consensus(d), f'{d}_consensus')


if __name__ == '__main__':
    main()
