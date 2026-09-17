#!/usr/bin/env python3
"""B04a_variant_density.py -- per-10kb-bin variant density from the assembly-vs-hg38 diploid VCFs
(G01, data/vcf_tmp/tmp_<chrom>/<sample>/hap_vs_hg38/diploid.vcf.gz), for testing whether PMD/domain
bins carry more sequence divergence than the rest of the genome.

Uses the per-chrom G01 output directly: the merged per-donor VCFs in data/vcf (G02/G03) were
built when only chr21 had finished, so they are chr21-only.

Counted per donor per bin, then summed over donors:
  snv, indel (<50bp), sv (>=50bp REF/ALT length difference), and het (GT 0|1 or 1|0) for each.
Also per-donor genome-wide totals, for normalisation.

Caveat kept in mind when interpreting: these are assembly-alignment calls, so satellite and
segmental-duplication sequence is under-called and repeat-rich bins are noisier. The companion
aggregation (QC15) controls for CpG count, gene content and mappability-like covariates.

Output: results/meth_bins/variant_density/<chrom>.bins.tsv.gz (bin_start, n_donors, snv, indel,
sv, snv_het, indel_het, sv_het) and <chrom>.per_donor.tsv.gz

Usage: B04a_variant_density.py <chrom>
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
VCF_TMP = PROJDIR / 'data' / 'vcf_tmp'
OUT = PROJDIR / 'results' / 'meth_bins' / 'variant_density'
CHROMSIZES = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosomal.chromsizes')
BIN = 10_000
SV_MIN = 50


def count_donor(path: Path, nbins: int) -> tuple:
    cols = {k: np.zeros(nbins, dtype=np.int32) for k in
            ('snv', 'indel', 'sv', 'snv_het', 'indel_het', 'sv_het')}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.split('\t', 10)
            pos = int(f[1])
            b = pos // BIN
            if b >= nbins:
                continue
            d = max(len(a) for a in f[4].split(',')) - len(f[3])
            kind = 'snv' if d == 0 and len(f[3]) == 1 else ('sv' if abs(d) >= SV_MIN else 'indel')
            cols[kind][b] += 1
            gt = f[9][:3]
            if gt in ('0|1', '1|0', '0/1', '1/0'):
                cols[kind + '_het'][b] += 1
    return cols


def main() -> None:
    chrom = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = pd.read_csv(CHROMSIZES, sep='\t', header=None, names=['chrom', 'len'])
    nbins = int(sizes.loc[sizes.chrom == chrom, 'len'].iloc[0] // BIN) + 1
    vcfs = sorted((VCF_TMP / f'tmp_{chrom}').glob('*/hap_vs_hg38/diploid.vcf.gz'))
    print(f'{chrom}: {len(vcfs)} donor VCFs, {nbins} bins', flush=True)
    total = {k: np.zeros(nbins, dtype=np.int64) for k in
             ('snv', 'indel', 'sv', 'snv_het', 'indel_het', 'sv_het')}
    per_donor = []
    for i, v in enumerate(vcfs):
        sample = v.parts[-3]
        c = count_donor(v, nbins)
        for k in total:
            total[k] += c[k]
        per_donor.append({'sample': sample, 'chrom': chrom, **{k: int(c[k].sum()) for k in c}})
        if (i + 1) % 25 == 0:
            print(f'  {i + 1}/{len(vcfs)}', flush=True)
    df = pd.DataFrame({'bin_start': np.arange(nbins) * BIN, 'n_donors': len(vcfs), **total})
    df.insert(0, 'chrom', chrom)
    df.to_csv(OUT / f'{chrom}.bins.tsv.gz', sep='\t', index=False)
    pd.DataFrame(per_donor).to_csv(OUT / f'{chrom}.per_donor.tsv.gz', sep='\t', index=False)
    print(f'# {chrom}: {df.snv.sum():,} SNV, {df.indel.sum():,} indel, {df.sv.sum():,} SV calls', flush=True)


if __name__ == '__main__':
    main()
