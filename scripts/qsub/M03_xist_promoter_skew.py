#!/usr/bin/env python3
"""M03_xist_promoter_skew.py -- per-haplotype methylation of the XIST promoter CpG island, as an
XCI-skew / clonality metric for female donors.

QC05 measured the whole ~32kb XIST gene body and found a muted, unimodal signal; the diagnostic
window is the promoter CpG island. XIST (ENSG00000229807, GRCh38 chrX:73,820,649-73,852,714,
minus strand) is silenced on the active X by promoter methylation and expressed from the inactive
X, so in a polyclonal population both haplotypes read ~50% methylated, while in a clonal
(skewed) population one haplotype is ~fully methylated and the other ~unmethylated.
skew = |hap1 mean - hap2 mean|.

The window is the hg38 CpG island overlapping the XIST TSS (cpgIslandExt), mapped to each
haplotype assembly with the existing chain files; reads come straight from the modbed, filtered
to positions that are CG in that assembly (same convention as A01a/C02a).

Output: results/qc/data/xist_promoter_skew.tsv
  sample, sex, hap1_meth, hap2_meth, skew, hap1_reads, hap2_reads, n_cpg_hap1, n_cpg_hap2

Usage: M03_xist_promoter_skew.py [n_workers]
"""
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'harmonize_hg38' / 'pilot'))
import chainmap  # noqa: E402

MODBED_DIR = PROJDIR / 'data' / 'modbed'
ASM_DIR = PROJDIR / 'data' / 'assemblies'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
COVS = PROJDIR / 'results' / 'qc' / 'data' / 'global_methylation_covariates_full_cohort.tsv'
CGI = Path('/u/project/cluo/terencew/reference/hg38/bed/cpgIslandExt.hg38.bed')
OUT = PROJDIR / 'results' / 'qc' / 'data' / 'xist_promoter_skew.tsv'
# XIST TSS from gencode v43 (chrX:73,820,649-73,852,723, minus strand), so the promoter lies at
# the high-coordinate end. The local cpgIslandExt track annotates no island here (its chrX entries
# around XIST look unreliable), so the window was placed empirically: CpG density measured in
# 500bp steps through an assembly peaks from TSS-1500 to TSS (~49 CpGs / 1.5kb) and drops to
# 2-5 per 500bp outside it.
CHROM, START, END = 'chrX', 73_851_123, 73_852_923


def hap_meth(sample: str, hap: int) -> tuple:
    cm = chainmap.ChainMap(str(CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'))
    tbx = pysam.TabixFile(str(MODBED_DIR / f'{sample}_hap{hap}.modbed.gz'))
    fa = pysam.FastaFile(str(ASM_DIR / f'{sample}_hap{hap}.fa.gz'))
    n_meth = n_tot = n_reads = 0
    positions = set()
    for contig, s, e in cm.q_interval_to_t(CHROM, START, END):
        full = f'{sample}#{hap}#{contig}'
        if full not in tbx.contigs:
            continue
        seq = np.frombuffer(fa.fetch(full, max(0, s - 2), e + 2).upper().encode(), dtype=np.uint8)
        off = max(0, s - 2)
        for row in tbx.fetch(full, s, e):
            f = row.split('\t')
            rs = int(f[1])
            shift = -1 if f[5] == '-' else 0
            n_reads += 1
            for fld, is_meth in ((f[6], 1), (f[7], 0)):
                if fld in ('', '.', '-'):
                    continue
                for o in fld.split(','):
                    p = rs + abs(int(o)) + shift
                    if not (s <= p < e):
                        continue
                    i = p - off
                    if not (0 <= i < len(seq) - 1 and seq[i] == ord('C') and seq[i + 1] == ord('G')):
                        continue
                    positions.add(p)
                    n_tot += 1
                    n_meth += is_meth
    return (n_meth / n_tot if n_tot else np.nan), n_reads, len(positions)


def run(sample: str) -> dict:
    try:
        m1, r1, c1 = hap_meth(sample, 1)
        m2, r2, c2 = hap_meth(sample, 2)
    except (FileNotFoundError, OSError, ValueError) as err:
        return {'sample': sample, 'error': str(err)[:80]}
    return {'sample': sample, 'hap1_meth': m1, 'hap2_meth': m2, 'skew': abs(m1 - m2),
            'hap1_reads': r1, 'hap2_reads': r2, 'n_cpg_hap1': c1, 'n_cpg_hap2': c2}


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    cov = pd.read_csv(COVS, sep='\t').groupby('sample_id').sex.first()
    females = [s for s, v in cov.items() if v == 'female'
               and (CHAIN_DIR / f'{s}_hap1_vs_GRCh38.chain.gz').exists()
               and (MODBED_DIR / f'{s}_hap1.modbed.gz').exists()]
    print(f'window {CHROM}:{START}-{END} ({END - START} bp); {len(females)} female donors', flush=True)
    with ProcessPoolExecutor(n) as ex:
        rows = list(ex.map(run, females))
    df = pd.DataFrame(rows)
    df['sex'] = 'female'
    df.to_csv(OUT, sep='\t', index=False)
    ok = df.dropna(subset=['skew']) if 'skew' in df else df
    if len(ok):
        print(ok[['hap1_meth', 'hap2_meth', 'skew']].describe().round(3).to_string())
        print(f"donors with skew > 0.5: {(ok['skew'] > 0.5).sum()} / {len(ok)}")


if __name__ == '__main__':
    main()
