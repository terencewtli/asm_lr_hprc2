#!/usr/bin/env python3
"""A02a_methcounts_to_bigwig.py -- per-haplotype hg38 CpG methylation bigWig for genome-browser
viewing (UCSC track hub: /u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2).

Data source: A01a's reference-CpG-filtered methcounts
(data/pmds/<sample>/<sample>_hap{N}.cpg.methcounts.tsv.gz), NOT P03's het-filtered matrix. The
k>=1 filter drops reads that don't span a het site, which would show filter-driven gaps as if
they were biology.

Coordinates: mapped to hg38 point-by-point with the vectorized chain lookup in
meth_bins/B01a_hap_bins.py (same convention, including the -1 shift on strand-flipped blocks;
validated against liftOver, 99.3% identical positions). The first version used UCSC liftOver on
~32M single-base records and every task hit the 1h h_rt limit (job 14772524, 2026-09-17).
hg38 CpGs hit twice are dropped; autosomes only (matches the chromsizes file).

Two bigWigs per haplotype: % methylation (0-100) and read depth.

Usage: A02a_methcounts_to_bigwig.py <sample> <hap> <out_meth.bw> <out_depth.bw>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'meth_bins'))
from B01a_hap_bins import load_blocks, map_points  # noqa: E402

PMD_DIR = PROJDIR / 'data' / 'pmds'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
HG38_CHROMSIZES = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosomal.chromsizes')


def main() -> None:
    sample, hap, out_meth, out_depth = sys.argv[1:5]
    methcounts_gz = PMD_DIR / sample / f'{sample}_hap{hap}.cpg.methcounts.tsv.gz'
    chain_gz = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'
    if not methcounts_gz.exists() or not chain_gz.exists():
        print(f'# SKIP {sample} hap{hap}: missing methcounts or chain', file=sys.stderr)
        sys.exit(0)

    m = pd.read_csv(methcounts_gz, sep='\t', header=None, usecols=[0, 1, 4, 5],
                    names=['contig', 'pos', 'frac', 'n'],
                    dtype={'contig': 'category', 'pos': np.int64, 'frac': np.float32, 'n': np.int32})
    contig = m.contig.values.astype(str)
    pos = m.pos.values
    blocks = load_blocks(chain_gz)
    qchrom = np.full(len(m), '', dtype=object)
    qpos = np.full(len(m), -1, dtype=np.int64)
    for c in np.unique(contig):
        bare = c.rsplit('#', 1)[-1]
        if bare not in blocks:
            continue
        idx = np.flatnonzero(contig == c)
        ok, qn, qp = map_points(blocks[bare], pos[idx])
        qchrom[idx[ok]] = qn[ok]
        qpos[idx[ok]] = qp[ok]

    sizes = pd.read_csv(HG38_CHROMSIZES, sep='\t', header=None, names=['chrom', 'len'])
    h = pd.DataFrame({'chrom': qchrom, 'pos': qpos, 'pct': m.frac.values * 100, 'n': m.n.values})
    h = h[h.chrom.isin(set(sizes.chrom))]
    h = h[~h.duplicated(['chrom', 'pos'], keep=False)]
    order = {c: i for i, c in enumerate(sizes.chrom)}
    h = h.assign(o=h.chrom.map(order)).sort_values(['o', 'pos'])

    for out, col in ((out_meth, 'pct'), (out_depth, 'n')):
        bw = pyBigWig.open(out + '.partial', 'w')
        bw.addHeader(list(zip(sizes.chrom, sizes.len.astype(int))))
        for c, g in h.groupby('o', sort=True):
            bw.addEntries(sizes.chrom.iloc[c], g.pos.values.astype(np.int64).tolist(), span=1,
                          values=g[col].values.astype(np.float64).tolist())
        bw.close()
        Path(out + '.partial').rename(out)
    print(f'# wrote {out_meth} and {out_depth} ({len(h):,} hg38 CpGs)', file=sys.stderr)


if __name__ == '__main__':
    main()
