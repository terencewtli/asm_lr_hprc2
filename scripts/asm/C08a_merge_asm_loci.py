#!/usr/bin/env python3
"""C08a_merge_asm_loci.py -- measure ASM locus SIZE, which C01a's tiling cannot.

C01a cuts the genome into <=1000 bp pieces (MAX_GAP 500, MAX_SPAN 1000, MIN_CPG 5), so the
"818 bp / 7 CpG ASM region" quoted in RESULTS is the tiling, not a property of ASM. This script
merges each donor's *adjacent significant* regions (bedtools-merge semantics, SLOP bp of slack to
absorb boundary jitter) and reports the distribution of merged locus extents -- the quantity
asm_lr got via `bedtools merge -d 250` and this project has never computed.

Reported per donor and pooled: n merged loci, how many tiles each absorbed, span in bp, total
CpGs, and the fraction of merged loci that are >1 tile (i.e. where the 1 kb cap was binding).
Restricted to the discovery tier by default so the size estimate is not a lambda artefact.

Outputs (results/asm/model/):
  asm_locus_size_per_donor.tsv   one row per donor
  asm_locus_size_pooled.tsv      span/tile-count histogram over the pooled discovery-tier loci
  asm_merged_loci.tsv.gz         every merged locus (sample, chrom, start, end, n_tiles, n_cpg)

Usage: C08a_merge_asm_loci.py [--tier discovery|replication|all] [--slop 250]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
GENOME = PROJDIR / 'results' / 'asm' / 'genome'
TIERS = PROJDIR / 'results' / 'asm' / 'replication' / 'donor_tiers.tsv'
OUT = PROJDIR / 'results' / 'asm' / 'model'


def merge_one(df: pd.DataFrame, slop: int) -> pd.DataFrame:
    # bedtools-merge over one donor's significant tiles: a new locus starts when the next tile
    # begins more than `slop` bp after the running end
    out = []
    for chrom, g in df.groupby('chrom'):
        g = g.sort_values('start')
        s, e, c = g.start.values, g.end.values, g.n_cpg_ref.values
        brk = np.r_[True, s[1:] > (np.maximum.accumulate(e)[:-1] + slop)]
        lid = np.cumsum(brk) - 1
        out.append(pd.DataFrame({
            'chrom': chrom,
            'start': np.bincount(lid, weights=s).astype(np.int64) // np.bincount(lid) * 0
                     + pd.Series(s).groupby(lid).min().values,
            'end': pd.Series(e).groupby(lid).max().values,
            'n_tiles': np.bincount(lid),
            'n_cpg': np.bincount(lid, weights=c).astype(np.int64)}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--tier', default='discovery')
    ap.add_argument('--slop', type=int, default=250)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    tiers = pd.read_csv(TIERS, sep='\t')
    samples = tiers['sample'] if a.tier == 'all' else tiers[tiers.tier == a.tier]['sample']
    print('tier=%s slop=%d -> %d donors' % (a.tier, a.slop, len(samples)), flush=True)

    rows, allloci = [], []
    for i, s in enumerate(samples, 1):
        f = GENOME / 'per_donor' / ('%s.asm_sig.tsv.gz' % s)
        if not f.exists():
            print('  skip %s (no asm_sig)' % s, file=sys.stderr)
            continue
        d = pd.read_csv(f, sep='\t', usecols=['region_id', 'set', 'chrom', 'start', 'end', 'n_cpg_ref'])
        d = d[d['set'] == 'cpg']
        if d.empty:
            continue
        mg = merge_one(d, a.slop)
        mg['span'] = mg.end - mg.start
        mg.insert(0, 'sample', s)
        allloci.append(mg)
        rows.append({'sample': s, 'n_tiles_sig': len(d), 'n_merged_loci': len(mg),
                     'tiles_per_locus': len(d) / len(mg),
                     'frac_multi_tile': float((mg.n_tiles > 1).mean()),
                     'span_median': float(mg.span.median()), 'span_p75': float(mg.span.quantile(.75)),
                     'span_p95': float(mg.span.quantile(.95)), 'span_max': int(mg.span.max()),
                     'cpg_median': float(mg.n_cpg.median())})
        if i % 20 == 0:
            print('  %d/%d' % (i, len(samples)), flush=True)

    per = pd.DataFrame(rows).merge(tiers[['sample', 'superpopulation', 'chemistry', 'lambda_gc']], on='sample')
    per.to_csv(OUT / 'asm_locus_size_per_donor.tsv', sep='\t', index=False, float_format='%.5g')

    loci = pd.concat(allloci, ignore_index=True)
    loci.to_csv(OUT / 'asm_merged_loci.tsv.gz', sep='\t', index=False)

    bins = [0, 1000, 2000, 3000, 5000, 10000, 20000, np.inf]
    labs = ['<=1kb (1 tile)', '1-2kb', '2-3kb', '3-5kb', '5-10kb', '10-20kb', '>20kb']
    hist = (pd.cut(loci.span, bins, labels=labs).value_counts(normalize=True).reindex(labs)
            .rename('frac').to_frame())
    hist['n'] = pd.cut(loci.span, bins, labels=labs).value_counts().reindex(labs)
    hist.index.name = 'span_bin'
    hist.to_csv(OUT / 'asm_locus_size_pooled.tsv', sep='\t', float_format='%.5g')

    print('\n=== merged ASM loci, tier=%s (%d donors) ===' % (a.tier, len(per)))
    print('median tiles per locus      %.2f' % per.tiles_per_locus.median())
    print('frac loci spanning >1 tile  %.3f' % per.frac_multi_tile.median())
    print('median merged span          %.0f bp   (tile cap is 1000 bp)' % per.span_median.median())
    print('pooled span distribution:\n' + hist.to_string())
    print('\npooled: n=%d loci, span median %.0f, p95 %.0f, max %d' % (
        len(loci), loci.span.median(), loci.span.quantile(.95), loci.span.max()))
    print('\nwrote', OUT / 'asm_locus_size_per_donor.tsv')


if __name__ == '__main__':
    main()
