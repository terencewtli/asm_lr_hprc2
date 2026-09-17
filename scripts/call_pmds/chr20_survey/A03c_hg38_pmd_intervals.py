#!/usr/bin/env python3
"""A03c_hg38_pmd_intervals.py -- hg38 PMD intervals from A03a's lifted per-CpG tables.

Rather than lifting PMD intervals (large intervals lift poorly and fragment), each run of
consecutive hg38 CpGs flagged in_own_pmd becomes one interval [first CpG, last CpG + 2); the run
breaks at any unflagged CpG. Runs shorter than MIN_CPG CpGs are dropped (single lifted CpGs
landing out of order). Done for every (sample, hap) with >=1 own PMD, plus a per-donor
consensus (CpGs flagged on both haps).

Outputs (results/qc/data/chr20_pmd_survey/hg38_pmds/):
  <sample>_hap<hap>.hg38.pmd.bed, <sample>_consensus.hg38.pmd.bed

Usage: A03c_hg38_pmd_intervals.py
"""
from pathlib import Path

import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
SURVEY = PROJDIR / 'results' / 'qc' / 'data' / 'chr20_pmd_survey'
OUT = SURVEY / 'hg38_pmds'
CHROM = 'chr20'
MIN_CPG = 20


def runs_to_bed(pos: pd.Series, flag: pd.Series) -> pd.DataFrame:
    df = pd.DataFrame({'pos': pos.values, 'flag': flag.values}).sort_values('pos')
    run_id = (df.flag != df.flag.shift()).cumsum()
    runs = df[df.flag == 1].groupby(run_id[df.flag == 1]).agg(
        start=('pos', 'min'), end=('pos', 'max'), n=('pos', 'size'))
    runs = runs[runs.n >= MIN_CPG]
    return pd.DataFrame({'chrom': CHROM, 'start': runs.start, 'end': runs.end + 2, 'n_cpg': runs.n})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    per_hap = pd.read_csv(SURVEY / 'per_hap.tsv', sep='\t')
    pos_samples = sorted(per_hap.loc[per_hap.n_pmd > 0, 'sample'].unique())
    for sample in pos_samples:
        tabs = {}
        for hap in (1, 2):
            f = SURVEY / 'per_hap' / f'{sample}_hap{hap}.hg38_cpg.tsv.gz'
            if not f.exists():
                continue
            t = pd.read_csv(f, sep='\t')
            tabs[hap] = t
            bed = runs_to_bed(t.pos, t.in_own_pmd)
            bed.to_csv(OUT / f'{sample}_hap{hap}.hg38.pmd.bed', sep='\t', header=False, index=False)
            print(f'{sample} hap{hap}: {len(bed)} intervals, {(bed.end - bed.start).sum() / 1e6:.1f} Mb')
        if len(tabs) == 2:
            m = tabs[1].merge(tabs[2], on='pos')
            flag = ((m.in_own_pmd_x == 1) & (m.in_own_pmd_y == 1)).astype(int)
            bed = runs_to_bed(m.pos, flag)
            bed.to_csv(OUT / f'{sample}_consensus.hg38.pmd.bed', sep='\t', header=False, index=False)
            print(f'{sample} consensus: {len(bed)} intervals, {(bed.end - bed.start).sum() / 1e6:.1f} Mb')


if __name__ == '__main__':
    main()
