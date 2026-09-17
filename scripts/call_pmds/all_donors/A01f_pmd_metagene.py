#!/usr/bin/env python3
"""A01f_pmd_metagene.py -- PMD metagene profiles for every donor haplotype over the shared hg38
reference sets from A01e.

Same idea as igvf 2023_YR2 pmds_by_rna_cluster_donor/03a-03c (PMD body split into scaled
windows plus fixed flanks, averaged over PMDs), but read directly from the per-CpG hg38 bigWigs
(A02a, % mCG) with pyBigWig.stats(nBins=...), so no window BED or tabix pass is needed.

Per PMD: BODY_BINS scaled bins over the body, FLANK_BINS fixed bins over FLANK bp each side.
Flank bins that overlap another interval of the same set are set to NaN (the old notebook
clipped flanks at the neighbouring PMD). Bin value = mean % mCG of CpGs in the bin
(unweighted); empty bins are NaN.

Output per (sample, hap): results/pmd_metagene/per_hap/<tag>.npz with, for each set,
  <set> -> float16 array (n_intervals, FLANK_BINS + BODY_BINS + FLANK_BINS)
Aggregation and plots: notebooks/qc/QC12_pmd_size_overlap_and_metagene.ipynb

Usage: A01f_pmd_metagene.py [n_workers]   (processes every *.hg38.meth.bw not yet done)
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
PMD_DIR = PROJDIR / 'data' / 'pmds'
REF = PROJDIR / 'results' / 'pmd_metagene' / 'ref_sets'
OUT = PROJDIR / 'results' / 'pmd_metagene' / 'per_hap'
SETS = ['lcl_constitutive_ge300kb', 'fibroblast', 'fibroblast_ge300kb', 'NA19338_consensus_ge300kb',
        'HG04187_consensus_ge300kb']
BODY_BINS = 40
FLANK_BINS = 20
FLANK = 200_000


def load_set(name: str, chrom_len: dict) -> pd.DataFrame:
    b = pd.read_csv(REF / f'{name}.bed', sep='\t', header=None, names=['chrom', 'start', 'end'])
    b = b[b.chrom.isin(chrom_len)].reset_index(drop=True)
    prev_end = b.groupby('chrom').end.shift(1).fillna(-1).values
    next_start = b.groupby('chrom').start.shift(-1).fillna(np.inf).values
    b['prev_end'], b['next_start'] = prev_end, next_start
    b['len'] = b.chrom.map(chrom_len)
    return b


def flank_mask(lo: int, hi: int, prev_end: float, next_start: float) -> np.ndarray:
    edges = np.linspace(lo, hi, FLANK_BINS + 1)
    return (edges[:-1] >= prev_end) & (edges[1:] <= next_start)


def profile(bw, r) -> np.ndarray:
    out = np.full(2 * FLANK_BINS + BODY_BINS, np.nan)
    body = bw.stats(r.chrom, int(r.start), int(r.end), type='mean', nBins=BODY_BINS)
    out[FLANK_BINS:FLANK_BINS + BODY_BINS] = [np.nan if v is None else v for v in body]
    for side, (lo, hi) in enumerate(((r.start - FLANK, r.start), (r.end, r.end + FLANK))):
        if lo < 0 or hi > r.len:
            continue
        v = np.array([np.nan if x is None else x for x in
                      bw.stats(r.chrom, int(lo), int(hi), type='mean', nBins=FLANK_BINS)])
        v[~flank_mask(lo, hi, r.prev_end, r.next_start)] = np.nan
        sl = slice(0, FLANK_BINS) if side == 0 else slice(FLANK_BINS + BODY_BINS, None)
        out[sl] = v
    return out


def run(bw_path: str) -> str:
    tag = Path(bw_path).name[:-len('.hg38.meth.bw')]
    out = OUT / f'{tag}.npz'
    if out.exists():
        return f'{tag} exists'
    bw = pyBigWig.open(bw_path)
    chrom_len = bw.chroms()
    res = {}
    for s in SETS:
        b = load_set(s, chrom_len)
        res[s] = np.vstack([profile(bw, r) for r in b.itertuples()]).astype(np.float16)
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + '.partial.npz')
    np.savez_compressed(tmp, **res)
    tmp.rename(out)
    return f'{tag} done'


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    OUT.mkdir(parents=True, exist_ok=True)
    todo = sorted(str(p) for p in PMD_DIR.glob('*/*.hg38.meth.bw'))
    with ProcessPoolExecutor(n) as ex:
        for msg in ex.map(run, todo):
            print(msg, flush=True)


if __name__ == '__main__':
    main()
