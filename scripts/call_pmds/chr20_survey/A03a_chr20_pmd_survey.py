#!/usr/bin/env python3
"""A03a_chr20_pmd_survey.py -- per-(sample, hap) chr20 PMD survey, feeding A03b's cross-donor
summary ("how methylated are NA19338's PMDs in every other donor?").

Steps, all restricted to the haplotype's chr20 contig (the assembly contig with the most
chain-aligned score to hg38 chr20):
  1. reference-CpG-filtered methcounts via A01a_modbed_to_methcounts.py (same code as the
     genome-wide run, contig-restricted)
  2. dnmtools pmd on that contig -> the haplotype's own PMD calls (native coordinates)
  3. per-CpG table lifted to hg38 (liftOver, same chain orientation as A02a: target=assembly,
     query=hg38), with an in_own_pmd flag set in native coordinates before lifting

Outputs (results/qc/data/chr20_pmd_survey/per_hap/):
  <sample>_hap<hap>.hg38_cpg.tsv.gz  pos, frac, n, in_own_pmd (hg38 chr20, 0-based CpG start)
  <sample>_hap<hap>.pmd.bed          own PMD calls, native contig coordinates
  <sample>_hap<hap>.summary.tsv      one-row summary

Usage: A03a_chr20_pmd_survey.py <sample> <hap> <tmpdir>
"""
import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
CHAIN_DIR = PROJDIR / 'data' / 'chains'
A01A = PROJDIR / 'scripts' / 'call_pmds' / 'all_donors' / 'A01a_modbed_to_methcounts.py'
OUTDIR = PROJDIR / 'results' / 'qc' / 'data' / 'chr20_pmd_survey' / 'per_hap'
DNMTOOLS = '/u/home/t/terencew/bin/dnmtools'
LIFTOVER = '/u/project/cluo/terencew/programs/kentsrc/utils/liftOver'


def chr20_contig(chain_gz: Path) -> str:
    score: dict = {}
    with gzip.open(chain_gz, 'rt') as fh:
        for line in fh:
            if line.startswith('chain'):
                f = line.split()
                if f[7] == 'chr20':
                    score[f[2]] = score.get(f[2], 0) + int(f[1])
    return max(score, key=score.get)


def main() -> None:
    sample, hap, tmpdir = sys.argv[1], sys.argv[2], Path(sys.argv[3])
    tmpdir.mkdir(parents=True, exist_ok=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tag = f'{sample}_hap{hap}'

    chain_gz = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'
    bare = chr20_contig(chain_gz)
    contig = f'{sample}#{hap}#{bare}'
    print(f'{tag}: chr20 contig {contig}', file=sys.stderr, flush=True)

    meth = tmpdir / f'{tag}.chr20.cpg.methcounts.tsv.gz'
    subprocess.run(['python3', str(A01A), sample, hap, str(meth), contig], check=True)

    pmd_bed = OUTDIR / f'{tag}.pmd.bed'
    subprocess.run([DNMTOOLS, 'pmd', '-o', str(pmd_bed), '-s', '1', '-i', '1000', str(meth)],
                   check=True)

    m = pd.read_csv(meth, sep='\t', header=None, usecols=[1, 4, 5], names=['pos', 'frac', 'n'])
    if pmd_bed.stat().st_size > 0:
        b = pd.read_csv(pmd_bed, sep='\t', header=None, usecols=[1, 2], names=['s', 'e'])
        b = b.sort_values('s').reset_index(drop=True)
    else:
        b = pd.DataFrame({'s': [], 'e': []}, dtype=int)
    if len(b):
        i = np.searchsorted(b.s.values, m.pos.values, side='right') - 1
        m['in_own_pmd'] = (i >= 0) & (m.pos.values < b.e.values[np.clip(i, 0, None)])
    else:
        m['in_own_pmd'] = False
    m['in_own_pmd'] = m.in_own_pmd.astype(int)

    native = tmpdir / f'{tag}.native.bed'
    lifted = tmpdir / f'{tag}.hg38.bed'
    unmapped = tmpdir / f'{tag}.unmapped'
    chain_plain = tmpdir / f'{tag}.chain'
    pd.DataFrame({'c': bare, 's': m.pos, 'e': m.pos + 1, 'frac': m.frac, 'n': m.n,
                  'pmd': m.in_own_pmd}).to_csv(native, sep='\t', header=False, index=False)
    with gzip.open(chain_gz, 'rt') as fh, open(chain_plain, 'w') as out:
        out.write(fh.read())
    subprocess.run([LIFTOVER, '-bedPlus=3', str(native), str(chain_plain), str(lifted),
                    str(unmapped)], check=True)
    h = pd.read_csv(lifted, sep='\t', header=None, names=['c', 's', 'e', 'frac', 'n', 'pmd'])
    h = h[h.c == 'chr20'].drop_duplicates('s', keep=False).sort_values('s')
    h[['s', 'frac', 'n', 'pmd']].rename(columns={'s': 'pos', 'pmd': 'in_own_pmd'}).to_csv(
        OUTDIR / f'{tag}.hg38_cpg.tsv.gz', sep='\t', index=False)

    contig_len = int(subprocess.run(
        f"zcat {chain_gz} | awk '$1==\"chain\" && $3==\"{bare}\" {{print $4; exit}}'",
        shell=True, check=True, capture_output=True, text=True).stdout.strip())
    L = (b.e - b.s) if len(b) else pd.Series([], dtype=float)
    summ = {
        'sample': sample, 'hap': int(hap), 'contig': bare, 'contig_len': contig_len,
        'n_cpg': len(m), 'mean_depth': m.n.mean(), 'median_depth': m.n.median(),
        'n_pmd': len(b), 'pmd_bp': int(L.sum()), 'pmd_frac_contig': L.sum() / contig_len,
        'pmd_mean_kb': L.mean() / 1e3 if len(b) else np.nan,
        'pmd_median_kb': L.median() / 1e3 if len(b) else np.nan,
        'meth_all': m.frac.mean(),
        'meth_in_own_pmd': m.frac[m.in_own_pmd == 1].mean() if len(b) else np.nan,
        'meth_out_own_pmd': m.frac[m.in_own_pmd == 0].mean(),
        'n_cpg_hg38': len(h), 'lift_rate': len(h) / len(m),
    }
    pd.DataFrame([summ]).to_csv(OUTDIR / f'{tag}.summary.tsv', sep='\t', index=False)
    print(pd.DataFrame([summ]).T.to_string(), file=sys.stderr)

    for f in (meth, native, lifted, unmapped, chain_plain):
        if f.exists():
            f.unlink()


if __name__ == '__main__':
    main()
