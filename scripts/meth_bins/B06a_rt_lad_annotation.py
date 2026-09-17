#!/usr/bin/env python3
"""B06a_rt_lad_annotation.py -- hg38 10kb-bin replication timing and LAD annotation, the direct
test of the PMD mechanism (late-replicating, lamina-associated domains lose methylation with
cell division; Zhou 2018, Endicott 2022).

Sources (hg19, lifted to hg38 with UCSC hg19ToHg38):
  RT  -- ENCODE/UW Repli-seq wavelet-smoothed signal, GM12878 (an LCL, so the matched cell type):
         wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig. High = early replicating.
  LAD -- UCSC laminB1Lads (Guelen 2008 LaminB1 DamID, Tig3 fibroblasts; the standard cLAD set,
         not LCL-derived -- noted as a caveat).

Method: score hg19 10kb bins (mean RT, LAD overlap bp), lift the bin intervals to hg38 with
liftOver, and assign to hg38 10kb bins by the lifted midpoint (bins whose lift is split or
missing are dropped).

Output: results/meth_bins/annotations/rt_lad_10kb.tsv.gz (chrom, bin_start, rt, lad_bp)

Usage: B06a_rt_lad_annotation.py
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
REF = Path('/u/project/cluo_scratch/terencew/claude/asm_lr_hprc2/ref')
OUT = PROJDIR / 'results' / 'meth_bins' / 'annotations'
RT_BW = REF / 'wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig'
LAD = REF / 'laminB1Lads.txt.gz'
CHAIN = REF / 'hg19ToHg38.over.chain.gz'
LIFTOVER = '/u/project/cluo/terencew/programs/kentsrc/utils/liftOver'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
BIN = 10_000


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bw = pyBigWig.open(str(RT_BW))
    lad = pd.read_csv(LAD, sep='\t', header=None, usecols=[1, 2, 3], names=['chrom', 'start', 'end'])
    lad = lad[lad.chrom.isin(AUTOSOMES)]

    rows = []
    for c in AUTOSOMES:
        L = bw.chroms().get(c)
        if L is None:
            continue
        n = L // BIN + 1
        rt = bw.stats(c, 0, n * BIN if n * BIN <= L else L, type='mean', nBins=n)
        lad_bp = np.zeros(n)
        for s, e in zip(lad[lad.chrom == c].start, lad[lad.chrom == c].end):
            for b in range(s // BIN, min(e // BIN + 1, n)):
                lad_bp[b] += max(0, min(e, (b + 1) * BIN) - max(s, b * BIN))
        rows.append(pd.DataFrame({'chrom': c, 'start': np.arange(n) * BIN,
                                  'rt': [np.nan if v is None else v for v in rt], 'lad_bp': lad_bp}))
    hg19 = pd.concat(rows, ignore_index=True)
    hg19['end'] = hg19.start + BIN
    hg19['id'] = np.arange(len(hg19))
    tmp_in, tmp_out, tmp_un = (OUT / f'.hg19_bins{s}' for s in ('.bed', '.hg38.bed', '.unmapped'))
    hg19[['chrom', 'start', 'end', 'id']].to_csv(tmp_in, sep='\t', header=False, index=False)
    subprocess.run([LIFTOVER, '-bedPlus=3', str(tmp_in), str(CHAIN), str(tmp_out), str(tmp_un)], check=True)
    lifted = pd.read_csv(tmp_out, sep='\t', header=None, names=['chrom', 'start', 'end', 'id'])
    lifted = lifted[lifted.chrom.isin(AUTOSOMES)]
    keep = lifted.id.value_counts()
    lifted = lifted[lifted.id.isin(keep[keep == 1].index)]     # dropped if the bin lifted in pieces
    lifted['mid'] = (lifted.start + lifted.end) // 2
    lifted['bin_start'] = (lifted.mid // BIN) * BIN
    m = lifted.merge(hg19[['id', 'rt', 'lad_bp']], on='id')
    out = m.groupby(['chrom', 'bin_start']).agg(rt=('rt', 'mean'), lad_bp=('lad_bp', 'mean')).reset_index()
    out.round(4).to_csv(OUT / 'rt_lad_10kb.tsv.gz', sep='\t', index=False)
    for f in (tmp_in, tmp_out, tmp_un):
        f.unlink()
    print(f'{len(hg19):,} hg19 bins -> {len(out):,} hg38 bins with RT/LAD; '
          f'RT non-missing {out.rt.notna().mean():.3f}, LAD-covered bins {(out.lad_bp > 5000).mean():.3f}')


if __name__ == '__main__':
    main()
