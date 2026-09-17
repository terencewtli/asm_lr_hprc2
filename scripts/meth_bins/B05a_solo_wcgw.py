#!/usr/bin/env python3
"""B05a_solo_wcgw.py -- solo-WCGW mitotic-clock methylation per haplotype.

Zhou et al. 2018 (Nat Genet, PMID 29610480) showed that methylation at "solo-WCGW" CpGs inside
PMDs is the cleanest readout of cumulative cell divisions: a CpG in a W-CG-W context (W = A/T)
with no other CpG within FLANK bp on either side. Those sites lose methylation fastest with
division, so they sharpen the domain-depth phenotype this project measures.

Step 1 (annotate, once): scan hg38 autosomes for [AT]CG[AT] with no other CpG within +/-FLANK,
write positions per chromosome to results/meth_bins/solo_wcgw/sites.npz.
Step 2 (per haplotype): read the per-CpG hg38 bigWig (A02a) in chunks and average methylation at
solo-WCGW sites vs all CpGs, split by LCL domain class (B01b domain_frequency_10kb).

Output: results/meth_bins/solo_wcgw/per_hap/<tag>.tsv (one row: class x site-set means/counts)

Usage: B05a_solo_wcgw.py annotate
       B05a_solo_wcgw.py run [n_workers]
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
import pysam

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
OUT = PROJDIR / 'results' / 'meth_bins' / 'solo_wcgw'
PER_HAP = OUT / 'per_hap'
PMD_DIR = PROJDIR / 'data' / 'pmds'
FREQ = PROJDIR / 'results' / 'meth_bins' / 'domain_frequency_10kb.tsv.gz'
HG38 = '/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
FLANK = 35
BIN = 10_000
CHUNK = 10_000_000
CLASSES = ['never', 'rare', 'variable', 'common', 'constitutive']


def annotate() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fa = pysam.FastaFile(HG38)
    res = {}
    for c in AUTOSOMES:
        a = np.frombuffer(fa.fetch(c).upper().encode(), dtype=np.uint8)
        cg = np.flatnonzero((a[:-1] == ord('C')) & (a[1:] == ord('G')))
        w_left = np.isin(a[np.clip(cg - 1, 0, None)], [ord('A'), ord('T')])
        w_right = np.isin(a[np.clip(cg + 2, 0, len(a) - 1)], [ord('A'), ord('T')])
        wcgw = w_left & w_right
        # solo: nearest other CpG further than FLANK on both sides
        d_prev = np.r_[np.inf, np.diff(cg)]
        d_next = np.r_[np.diff(cg), np.inf]
        solo = (d_prev > FLANK) & (d_next > FLANK)
        res[c] = cg[wcgw & solo].astype(np.int64)
        print(f'{c}: {len(cg):,} CpGs, {int(wcgw.sum()):,} WCGW, {len(res[c]):,} solo-WCGW', flush=True)
    np.savez_compressed(OUT / 'sites.npz', **res)
    print(f'total solo-WCGW: {sum(len(v) for v in res.values()):,}')


def class_lookup() -> dict:
    f = pd.read_csv(FREQ, sep='\t', usecols=['chrom', 'bin_start', 'class_pmd'])
    code = {c: i for i, c in enumerate(CLASSES)}
    out = {}
    for c, g in f.groupby('chrom'):
        arr = np.full(int(g.bin_start.max()) // BIN + 1, -1, dtype=np.int8)
        arr[g.bin_start.values // BIN] = [code[x] for x in g.class_pmd]
        out[c] = arr
    return out


def run_hap(bw_path: str) -> str:
    tag = Path(bw_path).name[:-len('.hg38.meth.bw')]
    out = PER_HAP / f'{tag}.tsv'
    if out.exists():
        return f'{tag} exists'
    sites = np.load(OUT / 'sites.npz')
    cls = class_lookup()
    bw = pyBigWig.open(bw_path)
    n_cls = len(CLASSES)
    acc = {k: np.zeros((n_cls, 2)) for k in ('solo', 'all')}   # [class, (sum, count)]
    for c in AUTOSOMES:
        if c not in bw.chroms():
            continue
        L = bw.chroms()[c]
        s_all = sites[c]
        carr = cls.get(c)
        for st in range(0, L, CHUNK):
            en = min(st + CHUNK, L)
            v = np.array(bw.values(c, st, en), dtype=np.float32)
            ok = np.flatnonzero(~np.isnan(v))
            if not len(ok):
                continue
            pos = ok + st
            b = pos // BIN
            valid = b < len(carr)
            pos, ok, b = pos[valid], ok[valid], b[valid]
            k = carr[b]
            vals = v[ok]
            keep = k >= 0
            np.add.at(acc['all'][:, 0], k[keep], vals[keep])
            np.add.at(acc['all'][:, 1], k[keep], 1)
            sel = s_all[(s_all >= st) & (s_all < en)]
            if len(sel):
                vs = v[sel - st]
                bs = sel // BIN
                good = (~np.isnan(vs)) & (bs < len(carr))
                ks = carr[bs[good]]
                vs = vs[good]
                keep2 = ks >= 0
                np.add.at(acc['solo'][:, 0], ks[keep2], vs[keep2])
                np.add.at(acc['solo'][:, 1], ks[keep2], 1)
    rows = []
    for site_set, a in acc.items():
        for i, c in enumerate(CLASSES):
            rows.append({'tag': tag, 'sample': tag.rsplit('_hap', 1)[0], 'hap': int(tag[-1]),
                         'site_set': site_set, 'class_pmd': c,
                         'n_cpg': int(a[i, 1]), 'mean_mcg': a[i, 0] / a[i, 1] if a[i, 1] else np.nan})
    PER_HAP.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    return f'{tag} done'


def main() -> None:
    if sys.argv[1] == 'annotate':
        annotate()
        return
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    todo = sorted(str(p) for p in PMD_DIR.glob('*/*.hg38.meth.bw'))
    with ProcessPoolExecutor(n) as ex:
        for msg in ex.map(run_hap, todo):
            print(msg, flush=True)


if __name__ == '__main__':
    main()
