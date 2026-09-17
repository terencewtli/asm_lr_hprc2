#!/usr/bin/env python3
"""A01g_downsample_pmd.py -- is per-haplotype coverage limiting for PMD calling?

The natural objection to calling PMDs per haplotype is "pool the two haplotypes to double
coverage". This tests the premise directly: thin one haplotype's methcounts to lower depth by
binomial subsampling of the read counts at each CpG (keeping the methylation rate unbiased), run
the identical `dnmtools pmd`, and compare the calls to the full-depth calls.

If calls are stable down to ~10-15x, the observed per-haplotype depth (cohort median 29.7x,
min 15.6x, none <15x) is already past the point where more reads change the answer, and pooling
would only cost haplotype resolution. Pooling is also not free: the two haplotypes are assembled
separately, so pooling requires projecting both into a common reference first.

Output: results/pmd_downsample/<sample>_hap<hap>.<depth>x.pmd.bed and a summary tsv with, for
each target depth, the PMD burden and the Jaccard vs the full-depth calls.

Usage: A01g_downsample_pmd.py <sample> <hap> [depths, default 20,15,10,5]
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
PMD_DIR = PROJDIR / 'data' / 'pmds'
OUT = PROJDIR / 'results' / 'pmd_downsample'
DNMTOOLS = '/u/home/t/terencew/bin/dnmtools'
BEDTOOLS = '/u/home/t/terencew/bin/bedtools'


def jaccard(a: Path, b: Path) -> float:
    for f in (a, b):
        if not f.exists() or f.stat().st_size == 0:
            return float('nan')
    srt = []
    for f in (a, b):
        p = f.with_suffix('.sorted.bed')
        p.write_text(subprocess.run([BEDTOOLS, 'sort', '-i', str(f)], check=True,
                                    capture_output=True, text=True).stdout)
        srt.append(p)
    out = subprocess.run([BEDTOOLS, 'jaccard', '-a', str(srt[0]), '-b', str(srt[1])], check=True,
                         capture_output=True, text=True).stdout
    for p in srt:
        p.unlink()
    return float(out.strip().split('\n')[1].split('\t')[2])


def main() -> None:
    sample, hap = sys.argv[1], sys.argv[2]
    depths = [int(x) for x in (sys.argv[3].split(',') if len(sys.argv) > 3 else ['20', '15', '10', '5'])]
    OUT.mkdir(parents=True, exist_ok=True)
    tag = f'{sample}_hap{hap}'
    full_mc = PMD_DIR / sample / f'{tag}.cpg.methcounts.tsv.gz'
    full_bed = PMD_DIR / sample / f'{tag}.pmd.bed'

    m = pd.read_csv(full_mc, sep='\t', header=None, names=['contig', 'pos', 'strand', 'ctx', 'frac', 'n'])
    obs = m.n.mean()
    print(f'{tag}: {len(m):,} CpGs, mean depth {obs:.1f}', flush=True)
    rng = np.random.default_rng(1)
    rows = []
    for d in depths:
        if d >= obs:
            continue
        keep = rng.binomial(m.n.values, d / obs)          # thin coverage, unbiased in rate
        meth = rng.binomial(keep, m.frac.values)          # redraw methylated count at that depth
        sub = pd.DataFrame({'contig': m.contig, 'pos': m.pos, 'strand': '+', 'ctx': 'CpG',
                            'frac': np.divide(meth, keep, out=np.zeros(len(keep)), where=keep > 0),
                            'n': keep})
        sub = sub[sub.n > 0]
        mc = OUT / f'{tag}.{d}x.methcounts.tsv.gz'
        sub.to_csv(mc, sep='\t', header=False, index=False, float_format='%.6f')
        bed = OUT / f'{tag}.{d}x.pmd.bed'
        subprocess.run([DNMTOOLS, 'pmd', '-o', str(bed), '-s', '1', '-i', '1000', str(mc)], check=True)
        b = pd.read_csv(bed, sep='\t', header=None, usecols=[1, 2]) if bed.stat().st_size else pd.DataFrame({1: [], 2: []})
        rows.append({'sample': sample, 'hap': int(hap), 'target_depth': d,
                     'actual_mean_depth': float(sub.n.mean()), 'n_cpg': len(sub), 'n_pmd': len(b),
                     'pmd_gb': float((b[2] - b[1]).sum() / 1e9) if len(b) else 0.0,
                     'jaccard_vs_full': jaccard(bed, full_bed)})
        mc.unlink()
        print(rows[-1], flush=True)

    full = pd.read_csv(full_bed, sep='\t', header=None, usecols=[1, 2])
    rows.append({'sample': sample, 'hap': int(hap), 'target_depth': round(obs), 'actual_mean_depth': obs,
                 'n_cpg': len(m), 'n_pmd': len(full), 'pmd_gb': float((full[2] - full[1]).sum() / 1e9),
                 'jaccard_vs_full': 1.0})
    pd.DataFrame(rows).to_csv(OUT / f'{tag}.downsample_summary.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == '__main__':
    main()
