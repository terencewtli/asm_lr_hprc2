#!/usr/bin/env python3
"""C01a_build_regions.py -- one hg38 candidate-region table shared by every donor, so ASM calls
are directly comparable across donors (penetrance = fraction of tested donors with ASM at the
same region).

Regions:
  cpg     -- hg38 reference CpGs (autosomes) clustered by distance: a new cluster starts after a
             gap > MAX_GAP bp (matches ont_asm_caller.region.cluster_cpgs' 500bp and DSS's
             smoothing span). Clusters longer than MAX_SPAN are cut into consecutive <=MAX_SPAN
             pieces. Clusters with < MIN_CPG CpGs are dropped.
  zink    -- Zink et al. 2018 parent-of-origin DMRs (GRCh38, asm_lr/data), the positive-control
             set; tested as whole regions.

Output: results/asm/regions/regions_hg38.tsv.gz
  region_id, set, chrom, start, end (0-based half-open, end = last CpG + 2), n_cpg_ref

Usage: C01a_build_regions.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
OUT = PROJDIR / 'results' / 'asm' / 'regions' / 'regions_hg38.tsv.gz'
HG38 = '/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa'
ZINK = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr/data/zink_2018_supp5_pofo_dmrs.csv')
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
MAX_GAP = 500
MAX_SPAN = 1000
MIN_CPG = 5


def cpg_positions(seq: str) -> np.ndarray:
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    return np.flatnonzero((a[:-1] == ord('C')) & (a[1:] == ord('G')))


def cluster(chrom: str, pos: np.ndarray) -> pd.DataFrame:
    brk = np.r_[True, np.diff(pos) > MAX_GAP]
    cid = np.cumsum(brk)
    first = pos[np.searchsorted(cid, cid)]          # first CpG of each position's cluster
    piece = (pos - first) // MAX_SPAN               # cut long clusters into <=MAX_SPAN pieces
    key = cid.astype(np.int64) * 100_000 + piece
    df = pd.DataFrame({'key': key, 'pos': pos}).groupby('key').pos.agg(['min', 'max', 'size'])
    df = df[df['size'] >= MIN_CPG]
    return pd.DataFrame({'set': 'cpg', 'chrom': chrom, 'start': df['min'].values,
                         'end': df['max'].values + 2, 'n_cpg_ref': df['size'].values})


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fa = pysam.FastaFile(HG38)
    parts = []
    zink = pd.read_csv(ZINK, skiprows=2, usecols=['Chrom', 'peakStart', 'peakStop'])
    zink = zink[zink.Chrom.isin(AUTOSOMES)]
    for chrom in AUTOSOMES:
        pos = cpg_positions(fa.fetch(chrom))
        c = cluster(chrom, pos)
        parts.append(c)
        z = zink[zink.Chrom == chrom]
        if len(z):
            n = [int(np.sum((pos >= s) & (pos < e))) for s, e in zip(z.peakStart, z.peakStop)]
            parts.append(pd.DataFrame({'set': 'zink', 'chrom': chrom, 'start': z.peakStart.values,
                                       'end': z.peakStop.values, 'n_cpg_ref': n}))
        print(f'{chrom}: {len(pos):,} CpGs -> {len(c):,} regions', flush=True)
    reg = pd.concat(parts, ignore_index=True)
    reg.insert(0, 'region_id', [f'{s}_{c}_{a}_{b}' for s, c, a, b in
                                zip(reg.set, reg.chrom, reg.start, reg.end)])
    reg.to_csv(OUT, sep='\t', index=False)
    print(reg.groupby('set').agg(n=('region_id', 'size'), cpgs=('n_cpg_ref', 'sum'),
                                 median_len=('end', lambda e: np.median(e - reg.loc[e.index, 'start']))))


if __name__ == '__main__':
    main()
