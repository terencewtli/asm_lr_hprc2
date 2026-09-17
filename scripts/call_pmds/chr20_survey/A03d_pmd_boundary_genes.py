#!/usr/bin/env python3
"""A03d_pmd_boundary_genes.py -- LCL (A03c hg38 PMDs) vs fibroblast (igvf_pgp) PMDs on chr20:
interval Jaccard, boundary concordance, and genes at boundaries.

- LCL intervals are first merged within MERGE_GAP (liftover can split one native PMD).
- Boundary = each PMD start/end; boundary window = +/- HALF_WIN (10kb total).
- A boundary is "shared" if its window overlaps a window of the other set (boundaries <=10kb
  apart), else set-specific.
- Genes (gencode v43 autosome dedup) overlapping boundary windows are listed per class; gene
  density per window is compared against SHUFFLES random 10kb windows on chr20.

Outputs (results/qc/data/chr20_pmd_survey/boundaries/):
  jaccard.tsv, boundaries_<lcl>.tsv, boundary_genes_<lcl>.tsv, lcl_boundary_recurrence.tsv

Usage: A03d_pmd_boundary_genes.py [lcl_sample ...]   (default: every non-empty consensus set)
"""
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
SURVEY = PROJDIR / 'results' / 'qc' / 'data' / 'chr20_pmd_survey'
PMD_DIR = SURVEY / 'hg38_pmds'
OUT = SURVEY / 'boundaries'
FIB = PROJDIR / 'reference' / 'igvf_pgp' / 'start_merged.filt_mcg.sorted.bed.gz'
GENES = Path('/u/project/cluo/terencew/reference/hg38_igvf/bed/genes/gencode.v43.autosome.dedup.names.bed.gz')
CHROM = 'chr20'
CHROM_LEN = 64444167
MERGE_GAP = 5000
HALF_WIN = 5000
SHUFFLES = 1000
BEDTOOLS = '/u/home/t/terencew/bin/bedtools'


def read_bed(path: Path, merge: bool) -> pd.DataFrame:
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
    df = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    df = df[df.chrom == CHROM].sort_values('start')
    if merge and len(df):
        grp = (df.start > df.end.cummax().shift() + MERGE_GAP).cumsum()
        df = df.groupby(grp).agg(chrom=('chrom', 'first'), start=('start', 'min'), end=('end', 'max'))
    return df.reset_index(drop=True)


def covered(df: pd.DataFrame) -> np.ndarray:
    mask = np.zeros(CHROM_LEN, dtype=bool)
    for s, e in zip(df.start, df.end):
        mask[s:e] = True
    return mask


def jaccard(a: pd.DataFrame, b: pd.DataFrame) -> tuple:
    ma, mb = covered(a), covered(b)
    inter, union = (ma & mb).sum(), (ma | mb).sum()
    return inter / union if union else np.nan, ma.sum(), mb.sum(), inter


def boundaries(df: pd.DataFrame, label: str) -> pd.DataFrame:
    b = pd.DataFrame({'pos': np.concatenate([df.start.values, df.end.values]),
                      'side': ['start'] * len(df) + ['end'] * len(df)})
    b['set'] = label
    return b.sort_values('pos').reset_index(drop=True)


def nearest(pos: np.ndarray, other: np.ndarray) -> np.ndarray:
    other = np.sort(other)
    i = np.clip(np.searchsorted(other, pos), 1, len(other) - 1)
    return np.minimum(np.abs(pos - other[i - 1]), np.abs(pos - other[i]))


def genes_in(windows: pd.DataFrame, genes: pd.DataFrame) -> list:
    out = []
    for s, e in zip(windows.pos - HALF_WIN, windows.pos + HALF_WIN):
        hit = genes[(genes.start < e) & (genes.end > s)]
        out.append(','.join(hit.name.tolist()))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fib = read_bed(FIB, merge=False)
    genes = pd.read_csv(GENES, sep='\t', header=None, usecols=[0, 1, 2, 3],
                        names=['chrom', 'start', 'end', 'name'])
    genes = genes[genes.chrom == CHROM]
    samples = sys.argv[1:] or sorted(
        p.name.split('_consensus')[0] for p in PMD_DIR.glob('*_consensus.hg38.pmd.bed')
        if p.stat().st_size > 0)
    lcl = {s: read_bed(PMD_DIR / f'{s}_consensus.hg38.pmd.bed', merge=True) for s in samples}

    rows = []
    for s, df in lcl.items():
        j, la, lb, i = jaccard(df, fib)
        rows.append({'a': s, 'b': 'fibroblast', 'jaccard': j, 'a_mb': la / 1e6, 'b_mb': lb / 1e6,
                     'inter_mb': i / 1e6, 'frac_a_in_b': i / la, 'frac_b_in_a': i / lb})
    for s1, s2 in combinations(lcl, 2):
        j, la, lb, i = jaccard(lcl[s1], lcl[s2])
        rows.append({'a': s1, 'b': s2, 'jaccard': j, 'a_mb': la / 1e6, 'b_mb': lb / 1e6,
                     'inter_mb': i / 1e6, 'frac_a_in_b': i / la, 'frac_b_in_a': i / lb})
    jac = pd.DataFrame(rows)
    jac.to_csv(OUT / 'jaccard.tsv', sep='\t', index=False)
    print(jac.round(3).to_string(index=False))

    fib_b = boundaries(fib, 'fibroblast')
    rng = np.random.default_rng(1)
    rand_pos = rng.integers(HALF_WIN, CHROM_LEN - HALF_WIN, size=SHUFFLES)
    rand_genes = np.array([len(g.split(',')) if g else 0
                           for g in genes_in(pd.DataFrame({'pos': rand_pos}), genes)])
    print(f'\nrandom 10kb windows: mean genes {rand_genes.mean():.2f}, '
          f'frac with >=1 gene {np.mean(rand_genes > 0):.2f}')

    all_lcl_b = []
    for s, df in lcl.items():
        b = boundaries(df, s)
        b['dist_to_fib'] = nearest(b.pos.values, fib_b.pos.values)
        b['class'] = np.where(b.dist_to_fib <= 2 * HALF_WIN, 'shared', 'lcl_specific')
        b['genes'] = genes_in(b, genes)
        b['n_genes'] = [len(g.split(',')) if g else 0 for g in b.genes]
        b.to_csv(OUT / f'boundaries_{s}.tsv', sep='\t', index=False)
        all_lcl_b.append(b)
        fb = fib_b.copy()
        fb['dist_to_lcl'] = nearest(fb.pos.values, b.pos.values)
        fb['class'] = np.where(fb.dist_to_lcl <= 2 * HALF_WIN, 'shared', 'fib_specific')
        fb['genes'] = genes_in(fb, genes)
        fb['n_genes'] = [len(g.split(',')) if g else 0 for g in fb.genes]
        both = pd.concat([b.drop(columns='dist_to_fib'), fb.drop(columns='dist_to_lcl')])
        both.to_csv(OUT / f'boundary_genes_{s}.tsv', sep='\t', index=False)
        print(f'\n{s}: {len(b)} LCL boundaries (median dist to nearest fibroblast boundary '
              f'{np.median(b.dist_to_fib) / 1e3:.0f} kb); {len(fb)} fibroblast boundaries')
        print(both.groupby('class').n_genes.agg(['size', 'mean', lambda x: np.mean(x > 0)])
              .rename(columns={'size': 'n_boundaries', 'mean': 'mean_genes', '<lambda_0>': 'frac_with_gene'})
              .round(2).to_string())

    # boundary recurrence across LCL donors: for each LCL boundary, how many donors have a
    # boundary within 10kb
    allb = pd.concat(all_lcl_b, ignore_index=True)
    rec = []
    for _, r in allb.iterrows():
        near = allb[(allb.set != r.set) & ((allb.pos - r.pos).abs() <= 2 * HALF_WIN)].set.nunique()
        rec.append(near + 1)
    allb['n_donors_with_boundary'] = rec
    allb.to_csv(OUT / 'lcl_boundary_recurrence.tsv', sep='\t', index=False)
    print(f'\nLCL boundary recurrence (donors with a boundary within 10kb, of {len(lcl)}):')
    print(allb.groupby('set').n_donors_with_boundary.describe()[['count', 'mean', '50%']].round(2).to_string())
    print('\nshared-vs-specific by recurrence (all LCL boundaries):')
    print(pd.crosstab(allb.n_donors_with_boundary, allb['class']).to_string())


if __name__ == '__main__':
    main()
