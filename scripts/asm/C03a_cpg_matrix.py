#!/usr/bin/env python3
"""C03a_cpg_matrix.py -- per-chromosome donor-haplotype x CpG matrices in hg38, over the UNION of
CpGs seen in any haplotype.

Two sources, written separately:
  hetfilt -- C02a cpg/<sample>_<chrom>.cpg.tsv.gz (reads spanning >=1 het site, HMMFlagger-
             filtered). This is the matrix matching the ASM calls.
  all     -- A02a bigWigs data/pmds/<s>/<s>_hap<h>.hg38.{meth,depth}.bw (all reads, from the
             CpG-filtered methcounts). n_meth = round(pct/100 * depth).

Output: results/asm/cpg_matrix/<chrom>.<source>.npz with
  pos (int64, sorted hg38 0-based CpG start), cols (str, '<sample>_hap<h>'),
  n_meth, n_total (uint16, shape n_cpg x n_cols; 0 total = not covered)

Usage: C03a_cpg_matrix.py <chrom> <hetfilt|all>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
CPG_DIR = PROJDIR / 'results' / 'asm' / 'cpg'
PMD_DIR = PROJDIR / 'data' / 'pmds'
OUT = PROJDIR / 'results' / 'asm' / 'cpg_matrix'


def load_hetfilt(chrom: str) -> dict:
    out = {}
    for f in sorted(CPG_DIR.glob(f'*_{chrom}.cpg.tsv.gz')):
        sample = f.name[:-len(f'_{chrom}.cpg.tsv.gz')]
        d = pd.read_csv(f, sep='\t')
        for hap, g in d.groupby('hap'):
            out[f'{sample}_hap{hap}'] = (g.pos.values, g.n_meth.values, (g.n_meth + g.n_unmeth).values)
    return out


def load_all(chrom: str) -> dict:
    import pyBigWig
    out = {}
    for f in sorted(PMD_DIR.glob('*/*.hg38.meth.bw')):
        tag = f.name[:-len('.hg38.meth.bw')]
        dpath = f.with_name(f'{tag}.hg38.depth.bw')
        if not dpath.exists():
            continue
        bm, bd = pyBigWig.open(str(f)), pyBigWig.open(str(dpath))
        if chrom not in bm.chroms():
            continue
        m = bm.intervals(chrom) or []
        dep = bd.intervals(chrom) or []
        pos = np.array([x[0] for x in m], dtype=np.int64)
        pct = np.array([x[2] for x in m], dtype=np.float64)
        tot = np.array([x[2] for x in dep], dtype=np.float64)
        out[tag] = (pos, np.rint(pct / 100 * tot), tot)
    return out


def main() -> None:
    chrom, source = sys.argv[1], sys.argv[2]
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_hetfilt(chrom) if source == 'hetfilt' else load_all(chrom)
    cols = sorted(data)
    pos = np.unique(np.concatenate([v[0] for v in data.values()]))
    n_meth = np.zeros((len(pos), len(cols)), dtype=np.uint16)
    n_tot = np.zeros((len(pos), len(cols)), dtype=np.uint16)
    for j, c in enumerate(cols):
        p, m, t = data.pop(c)
        i = np.searchsorted(pos, p)
        n_meth[i, j] = np.clip(m, 0, 65535)
        n_tot[i, j] = np.clip(t, 0, 65535)
    np.savez_compressed(OUT / f'{chrom}.{source}.npz', pos=pos, cols=np.array(cols),
                        n_meth=n_meth, n_total=n_tot)
    cov = (n_tot > 0).sum(axis=1)
    print(f'{chrom} {source}: {len(cols)} haps, {len(pos):,} union CpGs; covered in >=90% of haps: '
          f'{np.mean(cov >= 0.9 * len(cols)):.3f}; median depth (covered) '
          f'{np.median(n_tot[n_tot > 0]):.0f}', flush=True)


if __name__ == '__main__':
    main()
