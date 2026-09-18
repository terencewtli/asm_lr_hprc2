#!/usr/bin/env python3
"""C07a_asm_candidates.py -- ASM replication, step 1: donor tiers, candidate loci, and the static
(genotype-independent) annotation of every candidate.

Design (RESULTS §7 recipe): discover on well-calibrated donors, re-test everywhere else.
  - donor tiers from the chr20 genomic-control lambda (results/qc/data/asm_null_vs_real_chr20.tsv):
      discovery    lambda_gc <  1.2   (~69 donors; calls used as-is)
      replication  1.2 <= lambda <= 3 (~105 donors; genomic control applied in C07b)
      inflated     3 < lambda <= 8    (kept in outputs, excluded from penetrance)
      excluded     lambda > 8, or no lambda (NA19338, NA20806 have no null run)
  - candidates: every `cpg` region called ASM (C04a: BH q < 0.05 genome-wide, |delta| >= 0.2) in
    >= 1 discovery donor, plus all 229 Zink imprinted DMRs (set = zink, the positive control).

Static annotation per candidate (hg38):
  CpG architecture  width, n_cpg_ref, CpG per 100 bp, GC fraction, CpG obs/exp,
                    cgi_like (GC >= 0.5 and obs/exp >= 0.6, Gardiner-Garden on the region itself)
  imprinting        distance to the nearest Zink DMR; zink_10kb flag
  mQTL              overlaps an HPRC2 S11 promoter-mQTL CpG window (reference/hprc2)
  genes             distance to the nearest GENCODE v43 TSS + that gene's name/type
  chromatin         Vu & Ernst universal ChromHMM group with the largest overlap
  domains           10 kb bin PMD class / frequency (B01b), replication timing, LAD bp (B06a)
  discovery         n discovery donors with the call, mean |delta| among them

Outputs (results/asm/replication/):
  donor_tiers.tsv            sample, chemistry, superpopulation, lambda_gc, tier
  candidates.tsv.gz          one row per candidate, columns above

Usage: C07a_asm_candidates.py
"""
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import pysam

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
GENOME = PROJDIR / 'results' / 'asm' / 'genome'
REGIONS = PROJDIR / 'results' / 'asm' / 'regions' / 'regions_hg38.tsv.gz'
LAMBDA = PROJDIR / 'results' / 'qc' / 'data' / 'asm_null_vs_real_chr20.tsv'
MQTL = PROJDIR / 'reference' / 'hprc2' / 's11_promoter_mqtl_hg38.bed'
TSS = PROJDIR / 'reference' / 'hprc2' / 'gencode_v43_autosome_gene_tss.bed'
DOMAINS = PROJDIR / 'results' / 'meth_bins' / 'domain_frequency_10kb.tsv.gz'
RT_LAD = PROJDIR / 'results' / 'meth_bins' / 'annotations' / 'rt_lad_10kb.tsv.gz'
CHROMHMM = '/u/project/cluo/terencew/reference/hg38_igvf/bed/chromhmm/vu.states.bed.gz'
HG38 = '/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa'
OUT = PROJDIR / 'results' / 'asm' / 'replication'
BIN = 10_000


def donor_tiers() -> pd.DataFrame:
    ds = pd.read_csv(GENOME / 'donor_summary.tsv', sep='\t')[['sample', 'chemistry', 'superpopulation']]
    lam = pd.read_csv(LAMBDA, sep='\t')[['sample', 'lambda_gc']]
    d = ds.merge(lam, on='sample', how='left')
    d['tier'] = np.select(
        [d.lambda_gc.isna(), d.lambda_gc < 1.2, d.lambda_gc <= 3, d.lambda_gc <= 8],
        ['excluded', 'discovery', 'replication', 'inflated'], default='excluded')
    return d


def nearest(chrom_pos: pd.DataFrame, ref: pd.DataFrame, cols: Dict[str, str]) -> pd.DataFrame:
    # distance from each region [start, end) to the nearest ref interval [start, end) on the same
    # chrom (0 if overlapping); carries the named ref columns of that nearest interval
    out = pd.DataFrame(index=chrom_pos.index, columns=['dist'] + list(cols.values()), dtype=object)
    for chrom, g in chrom_pos.groupby('chrom'):
        r = ref[ref.chrom == chrom].sort_values('start')
        if r.empty:
            continue
        rs, re_ = r.start.values, r.end.values
        mid = ((g.start.values + g.end.values) // 2)
        i = np.searchsorted(rs, mid)
        best = np.zeros(len(g), dtype=np.int64)
        dist = np.full(len(g), np.inf)
        # neighbours on both sides of the midpoint (a long interval starting 2 back can still
        # overlap, so check i-2 too)
        for j in (np.clip(i - 2, 0, len(r) - 1), np.clip(i - 1, 0, len(r) - 1), np.clip(i, 0, len(r) - 1)):
            d = np.maximum(0, np.maximum(rs[j] - g.end.values, g.start.values - re_[j]))
            better = d < dist
            dist[better], best[better] = d[better], j[better]
        out.loc[g.index, 'dist'] = dist
        for src, dst in cols.items():
            out.loc[g.index, dst] = r[src].values[best]
    out['dist'] = out['dist'].astype(float)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tiers = donor_tiers()
    tiers.to_csv(OUT / 'donor_tiers.tsv', sep='\t', index=False)
    print('donor tiers:', tiers.tier.value_counts().to_dict(), flush=True)

    disc = tiers[tiers.tier == 'discovery']['sample']
    sig = pd.concat([pd.read_csv(GENOME / 'per_donor' / f'{s}.asm_sig.tsv.gz', sep='\t',
                                 usecols=['sample', 'region_id', 'set', 'delta']) for s in disc])
    sig = sig[sig.set == 'cpg']
    rec = sig.groupby('region_id').agg(n_asm_discovery=('sample', 'size'),
                                       mean_absdelta_discovery=('delta', lambda x: x.abs().mean()))
    reg = pd.read_csv(REGIONS, sep='\t')
    cand = reg[reg.region_id.isin(rec.index) | (reg.set == 'zink')].copy()
    cand = cand.merge(rec, left_on='region_id', right_index=True, how='left')
    cand['n_asm_discovery'] = cand.n_asm_discovery.fillna(0).astype(int)
    cand['n_discovery_donors'] = len(disc)
    print(f'{len(cand):,} candidates ({(cand.set == "cpg").sum():,} cpg from {len(disc)} discovery donors, '
          f'{(cand.set == "zink").sum()} zink)', flush=True)

    # CpG architecture from the reference sequence
    fa = pysam.FastaFile(HG38)
    gc, cpg, oe = [], [], []
    for c, s, e in zip(cand.chrom, cand.start, cand.end):
        seq = fa.fetch(c, int(s), int(e)).upper()
        n = max(len(seq), 1)
        nc, ng, ncg = seq.count('C'), seq.count('G'), seq.count('CG')
        gc.append((nc + ng) / n)
        cpg.append(100 * ncg / n)
        oe.append(ncg * n / (nc * ng) if nc and ng else 0.0)
    cand['width'] = cand.end - cand.start
    cand['gc_frac'], cand['cpg_per_100bp'], cand['cpg_obs_exp'] = gc, cpg, oe
    cand['cgi_like'] = (cand.gc_frac >= 0.5) & (cand.cpg_obs_exp >= 0.6)

    # imprinting: nearest Zink DMR
    zink = reg[reg.set == 'zink'][['chrom', 'start', 'end', 'region_id']].rename(columns={'region_id': 'zink_id'})
    z = nearest(cand[['chrom', 'start', 'end']], zink, {'zink_id': 'nearest_zink'})
    cand['dist_zink'], cand['nearest_zink'] = z.dist.values, z.nearest_zink.values
    cand['zink_10kb'] = cand.dist_zink <= 10_000

    # mQTL: HPRC2 S11 promoter mQTL windows (hg38 cols 1-3)
    mq = pd.read_csv(MQTL, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    m = nearest(cand[['chrom', 'start', 'end']], mq, {})
    cand['dist_mqtl'] = m.dist.values
    cand['mqtl_overlap'] = cand.dist_mqtl == 0

    # genes: nearest TSS
    tss = pd.read_csv(TSS, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene', 'gene_type', 'strand'])
    t = nearest(cand[['chrom', 'start', 'end']], tss, {'gene': 'nearest_gene', 'gene_type': 'nearest_gene_type'})
    cand['dist_tss'], cand['nearest_gene'], cand['nearest_gene_type'] = t.dist.values, t.nearest_gene.values, t.nearest_gene_type.values

    # chromatin: ChromHMM group with the largest overlap
    tb = pysam.TabixFile(CHROMHMM)
    grp = []
    for c, s, e in zip(cand.chrom, cand.start, cand.end):
        ov: Dict[str, int] = {}
        for line in tb.fetch(c, int(s), int(e)):
            f = line.split('\t')
            ov[f[6]] = ov.get(f[6], 0) + min(int(f[2]), e) - max(int(f[1]), s)
        grp.append(max(ov, key=ov.get) if ov else 'NA')
    cand['chromhmm_group'] = grp

    # domains / RT / LAD of the 10 kb bin holding the region midpoint
    cand['bin_start'] = ((cand.start + cand.end) // 2) // BIN * BIN
    dom = pd.read_csv(DOMAINS, sep='\t', usecols=['chrom', 'bin_start', 'freq_pmd', 'class_pmd'])
    rt = pd.read_csv(RT_LAD, sep='\t')
    cand = cand.merge(dom, on=['chrom', 'bin_start'], how='left').merge(rt, on=['chrom', 'bin_start'], how='left')

    cand.to_csv(OUT / 'candidates.tsv.gz', sep='\t', index=False, float_format='%.5g')
    c = cand[cand.set == 'cpg']
    print('\ncpg candidates by discovery recurrence:')
    print(pd.cut(c.n_asm_discovery, [0, 1, 5, 20, 34, 62, 69]).value_counts().sort_index().to_string())
    print(f'\nzink within 10kb: {c.zink_10kb.mean():.3f}; mQTL overlap: {c.mqtl_overlap.mean():.3f}; '
          f'cgi_like: {c.cgi_like.mean():.3f}; median dist to TSS: {c.dist_tss.median():.0f} bp')
    print('\nwrote', OUT / 'candidates.tsv.gz')


if __name__ == '__main__':
    main()
