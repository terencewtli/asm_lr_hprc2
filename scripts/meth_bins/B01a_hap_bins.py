#!/usr/bin/env python3
"""B01a_hap_bins.py -- per-(sample, hap) genome-wide hg38 10kb methylation bins + hg38 PMD
intervals, the input to B01b's sample x bin matrix / PCA / domain-frequency analyses.

Inputs (all whole-genome):
  data/pmds/<s>/<s>_hap<h>.cpg.methcounts.tsv.gz   A01a, reference-CpG-filtered, native coords
  data/pmds/<s>/<s>_hap<h>.pmd.bed                 A01c genome-wide dnmtools pmd (native)
  data/chains/<s>_hap<h>_vs_GRCh38.chain.gz        target=assembly, query=hg38

Each CpG is flagged in_pmd in native coordinates, then mapped to hg38 point-by-point with a
vectorized chain lookup (numpy searchsorted over ungapped blocks from
harmonize_hg38/pilot/chainmap.parse_chain). UCSC liftOver was ~3h per haplotype on ~32M
single-base records (measured 2026-09-17), too slow for 404 haps. Convention matches
ChainMap.t_to_q: the block with the largest start <= pos is used. Both bases of the CpG must fall
in the block. On strand-flipped blocks the native C maps to the hg38 G, so the hg38 CpG
coordinate is (mapped - 1). hg38 CpGs hit twice are dropped. Autosomes only.

Outputs (results/meth_bins/per_hap/):
  <tag>.bins10kb.tsv.gz  chrom, bin_start, n_cpg, mean_meth (per-CpG mean), wmean_meth
                         (depth-weighted), mean_depth, frac_pmd
  <tag>.hg38.pmd.bed     runs of consecutive in_pmd hg38 CpGs (>=MIN_RUN_CPG), merged within
                         MERGE_GAP
  <tag>.summary.tsv      native global meth / depth / PMD burden, mapping rate

Usage: B01a_hap_bins.py <sample> <hap>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'harmonize_hg38' / 'pilot'))
from chainmap import parse_chain  # noqa: E402

PMD_DIR = PROJDIR / 'data' / 'pmds'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
OUTDIR = PROJDIR / 'results' / 'meth_bins' / 'per_hap'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
BIN = 10_000
MIN_RUN_CPG = 20
MERGE_GAP = 5000


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def load_blocks(chain_gz: Path) -> dict:
    b = pd.DataFrame(parse_chain(chain_gz), columns=['t', 't0', 't1', 'q', 'q0', 'q1', 'flip'])
    return {t: g.sort_values('t0').reset_index(drop=True) for t, g in b.groupby('t')}


def map_points(blocks: pd.DataFrame, pos: np.ndarray) -> tuple:
    i = np.searchsorted(blocks.t0.values, pos, side='right') - 1
    ic = np.clip(i, 0, None)
    t0, t1 = blocks.t0.values[ic], blocks.t1.values[ic]
    ok = (i >= 0) & (pos >= t0) & (pos + 1 < t1)
    off = pos - t0
    flip = blocks.flip.values[ic]
    q = np.where(flip, blocks.q1.values[ic] - 1 - off - 1, blocks.q0.values[ic] + off)
    return ok, blocks.q.values[ic], q


def flag_pmd(contig: np.ndarray, pos: np.ndarray, pmd_bed: Path) -> np.ndarray:
    flag = np.zeros(len(pos), dtype=np.int8)
    if not pmd_bed.exists() or pmd_bed.stat().st_size == 0:
        return flag
    b = pd.read_csv(pmd_bed, sep='\t', header=None, usecols=[0, 1, 2], names=['contig', 's', 'e'])
    for c, bc in b.groupby('contig'):
        idx = np.flatnonzero(contig == c)
        if not len(idx):
            continue
        bc = bc.sort_values('s')
        p = pos[idx]
        i = np.searchsorted(bc.s.values, p, side='right') - 1
        inside = (i >= 0) & (p < bc.e.values[np.clip(i, 0, None)])
        flag[idx[inside]] = 1
    return flag


def pmd_runs(h: pd.DataFrame) -> pd.DataFrame:
    out = []
    for chrom, g in h.groupby('chrom', sort=False):
        f = g.pmd.values
        run_id = np.cumsum(np.r_[1, f[1:] != f[:-1]])
        sel = f == 1
        if not sel.any():
            continue
        r = pd.DataFrame({'run': run_id[sel], 'pos': g.pos.values[sel]}).groupby('run').pos.agg(
            ['min', 'max', 'size'])
        r = r[r['size'] >= MIN_RUN_CPG].sort_values('min')
        if not len(r):
            continue
        s, e = r['min'].values, r['max'].values + 2
        grp = np.cumsum(np.r_[1, s[1:] > np.maximum.accumulate(e)[:-1] + MERGE_GAP])
        mg = pd.DataFrame({'s': s, 'e': e, 'g': grp}).groupby('g').agg(s=('s', 'min'), e=('e', 'max'))
        out.append(pd.DataFrame({'chrom': chrom, 'start': mg.s.values, 'end': mg.e.values}))
    return pd.concat(out) if out else pd.DataFrame(columns=['chrom', 'start', 'end'])


def main() -> None:
    sample, hap = sys.argv[1], sys.argv[2]
    tag = f'{sample}_hap{hap}'
    OUTDIR.mkdir(parents=True, exist_ok=True)

    meth = PMD_DIR / sample / f'{tag}.cpg.methcounts.tsv.gz'
    pmd_bed = PMD_DIR / sample / f'{tag}.pmd.bed'
    chain_gz = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'

    log(f'{tag}: reading methcounts')
    m = pd.read_csv(meth, sep='\t', header=None, usecols=[0, 1, 4, 5],
                    names=['contig', 'pos', 'frac', 'n'],
                    dtype={'contig': 'category', 'pos': np.int64, 'frac': np.float32, 'n': np.int32})
    contig = m.contig.values.astype(str)
    pos = m.pos.values
    m['pmd'] = flag_pmd(contig, pos, pmd_bed)
    pmd_native_bp = 0
    if pmd_bed.exists() and pmd_bed.stat().st_size > 0:
        b = pd.read_csv(pmd_bed, sep='\t', header=None, usecols=[1, 2])
        pmd_native_bp = int((b[2] - b[1]).sum())
    summ = {'sample': sample, 'hap': int(hap), 'n_cpg_native': len(m),
            'global_meth': float(m.frac.mean()), 'global_wmeth': float((m.frac * m.n).sum() / m.n.sum()),
            'mean_depth': float(m.n.mean()), 'median_depth': float(m.n.median()),
            'pmd_native_bp': pmd_native_bp, 'frac_cpg_in_pmd_native': float(m.pmd.mean()),
            'meth_in_pmd_native': float(m.frac[m.pmd == 1].mean()) if m.pmd.any() else np.nan,
            'meth_out_pmd_native': float(m.frac[m.pmd == 0].mean())}

    log(f'{tag}: mapping {len(m):,} CpGs to hg38')
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
    h = pd.DataFrame({'chrom': qchrom, 'pos': qpos, 'frac': m.frac.values, 'n': m.n.values,
                      'pmd': m.pmd.values})
    del m
    h = h[h.chrom.isin(AUTOSOMES)]
    h = h.sort_values(['chrom', 'pos'])
    h = h[~h.duplicated(['chrom', 'pos'], keep=False)].reset_index(drop=True)
    summ['n_cpg_hg38'] = len(h)
    summ['map_rate'] = len(h) / summ['n_cpg_native']

    log(f'{tag}: binning')
    h['bin_start'] = (h.pos // BIN) * BIN
    h['fn'] = h.frac * h.n
    bins = h.groupby(['chrom', 'bin_start'], sort=False).agg(
        n_cpg=('pos', 'size'), mean_meth=('frac', 'mean'), fn=('fn', 'sum'), nsum=('n', 'sum'),
        mean_depth=('n', 'mean'), frac_pmd=('pmd', 'mean')).reset_index()
    bins['wmean_meth'] = bins.fn / bins.nsum
    bins = bins[['chrom', 'bin_start', 'n_cpg', 'mean_meth', 'wmean_meth', 'mean_depth', 'frac_pmd']]
    bins.round(4).to_csv(OUTDIR / f'{tag}.bins10kb.tsv.gz', sep='\t', index=False)

    runs = pmd_runs(h)
    runs.to_csv(OUTDIR / f'{tag}.hg38.pmd.bed', sep='\t', header=False, index=False)
    summ['pmd_hg38_bp'] = int((runs.end - runs.start).sum()) if len(runs) else 0
    summ['n_pmd_hg38'] = len(runs)
    pd.DataFrame([summ]).to_csv(OUTDIR / f'{tag}.summary.tsv', sep='\t', index=False)
    log(pd.DataFrame([summ]).T.to_string())


if __name__ == '__main__':
    main()
