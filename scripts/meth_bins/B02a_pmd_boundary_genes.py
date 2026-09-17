#!/usr/bin/env python3
"""B02a_pmd_boundary_genes.py -- genome-wide LCL vs fibroblast PMD boundaries and the genes at
them. Genome-wide successor of chr20_survey/A03d.

LCL PMDs per donor = intersection of hap1 and hap2 hg38 PMD intervals (B01a), merged within
MERGE_GAP. Boundary = interval start/end.
- Boundary frequency map: for every 10kb hg38 bin, the fraction of donors with a boundary in
  that bin or an adjacent one (+/-1 bin).
- Recurrent LCL boundary bins: frequency >= REC_FREQ, collapsed into loci (adjacent bins).
- Classes: shared = recurrent LCL locus within SHARE_DIST of a fibroblast boundary;
  lcl_specific = recurrent LCL locus with no fibroblast boundary nearby; fib_specific =
  fibroblast boundary with no recurrent LCL locus nearby.
- Genes (gencode v43 autosome dedup) overlapping +/-HALF_WIN of each boundary/locus center.
  Enrichment: fraction of windows with >=1 gene vs N_PERM sets of random windows placed
  per chromosome in the same numbers.

Outputs (results/meth_bins/boundaries/): donor_fib_jaccard.tsv, boundary_freq_10kb.tsv.gz,
boundary_loci.tsv, boundary_gene_enrichment.tsv, boundary_genes_by_class.tsv

Usage: B02a_pmd_boundary_genes.py
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
PER_HAP = PROJDIR / 'results' / 'meth_bins' / 'per_hap'
OUT = PROJDIR / 'results' / 'meth_bins' / 'boundaries'
FIB = PROJDIR / 'reference' / 'igvf_pgp' / 'start_merged.filt_mcg.sorted.bed.gz'
GENES = Path('/u/project/cluo/terencew/reference/hg38_igvf/bed/genes/gencode.v43.autosome.dedup.names.bed.gz')
CHROMSIZES = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosomal.chromsizes')
BEDTOOLS = '/u/home/t/terencew/bin/bedtools'
BIN = 10_000
MERGE_GAP = 5000
HALF_WIN = 5000
SHARE_DIST = 10_000
REC_FREQ = 0.5
N_PERM = 200


def log(msg: str) -> None:
    print(msg, flush=True)


def bt(args: list) -> str:
    return subprocess.run([BEDTOOLS] + args, check=True, capture_output=True, text=True).stdout


def read_bed_str(s: str) -> pd.DataFrame:
    rows = [l.split('\t')[:3] for l in s.strip().split('\n') if l]
    df = pd.DataFrame(rows, columns=['chrom', 'start', 'end'])
    return df.astype({'start': int, 'end': int})


def gene_hits(chrom: np.ndarray, center: np.ndarray, genes: dict) -> list:
    out = []
    for c, p in zip(chrom, center):
        g = genes.get(c)
        if g is None:
            out.append('')
            continue
        hit = g[(g.start < p + HALF_WIN) & (g.end > p - HALF_WIN)]
        out.append(','.join(hit.name))
    return out


def frac_with_gene(chrom: np.ndarray, center: np.ndarray, gstart: dict, gend_max: dict) -> float:
    n = 0
    for c in np.unique(chrom):
        p = center[chrom == c]
        s, emax = gstart[c], gend_max[c]
        # any gene with start < p+W and end > p-W: among genes with start < p+W, max end > p-W
        i = np.searchsorted(s, p + HALF_WIN) - 1
        n += np.sum((i >= 0) & (emax[np.clip(i, 0, None)] > p - HALF_WIN))
    return n / len(center)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = OUT / 'tmp'
    tmp.mkdir(exist_ok=True)
    sizes = pd.read_csv(CHROMSIZES, sep='\t', header=None, names=['chrom', 'len'])
    sizes = sizes[sizes.chrom.isin([f'chr{i}' for i in range(1, 23)])]
    nbins = dict(zip(sizes.chrom, (sizes.len + BIN - 1) // BIN))
    fib = pd.read_csv(FIB, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
    fib = fib[fib.chrom.isin(nbins)]
    fib_path = tmp / 'fib.bed'
    fib.sort_values(['chrom', 'start']).to_csv(fib_path, sep='\t', header=False, index=False)

    samples = sorted({f.name.split('_hap')[0] for f in PER_HAP.glob('*_hap1.hg38.pmd.bed')}
                     & {f.name.split('_hap')[0] for f in PER_HAP.glob('*_hap2.hg38.pmd.bed')})
    log(f'{len(samples)} donors with both haps')
    counts = {c: np.zeros(n, dtype=np.int32) for c, n in nbins.items()}
    jrows = []
    for s in samples:
        a, b = (PER_HAP / f'{s}_hap{h}.hg38.pmd.bed' for h in (1, 2))
        if a.stat().st_size == 0 or b.stat().st_size == 0:
            continue
        inter = bt(['intersect', '-a', str(a), '-b', str(b)])
        cons = tmp / f'{s}.bed'
        cons.write_text(inter)
        inter = bt(['sort', '-i', str(cons)]) if inter else ''
        cons.write_text(inter)
        merged = bt(['merge', '-d', str(MERGE_GAP), '-i', str(cons)]) if inter else ''
        cons.write_text(merged)
        df = read_bed_str(merged)
        df = df[df.chrom.isin(nbins)]
        j = bt(['jaccard', '-a', str(cons), '-b', str(fib_path)]).strip().split('\n')[1].split('\t')
        jrows.append({'sample': s, 'pmd_mb': (df.end - df.start).sum() / 1e6, 'n_pmd': len(df),
                      'jaccard_fib': float(j[2]), 'inter_mb': int(j[0]) / 1e6})
        for c, g in df.groupby('chrom'):
            bb = np.unique(np.concatenate([g.start.values, g.end.values]) // BIN)
            hit = np.unique(np.clip(np.concatenate([bb - 1, bb, bb + 1]), 0, nbins[c] - 1))
            counts[c][hit] += 1
        cons.unlink()
    jac = pd.DataFrame(jrows)
    jac.to_csv(OUT / 'donor_fib_jaccard.tsv', sep='\t', index=False)
    n = len(jac)
    log(f'genome-wide consensus PMD vs fibroblast Jaccard (n={n}): median {jac.jaccard_fib.median():.3f} '
        f'range {jac.jaccard_fib.min():.3f}-{jac.jaccard_fib.max():.3f}; PMD Mb median {jac.pmd_mb.median():.0f}')

    freq = pd.concat([pd.DataFrame({'chrom': c, 'bin_start': np.arange(len(v)) * BIN, 'freq': v / n})
                      for c, v in counts.items()], ignore_index=True)
    freq.to_csv(OUT / 'boundary_freq_10kb.tsv.gz', sep='\t', index=False)
    log('boundary-bin frequency quantiles: ' + freq.freq.quantile([.5, .9, .99, .999]).round(3).to_dict().__repr__())

    rec = freq[freq.freq >= REC_FREQ].copy()
    rec['locus'] = ((rec.chrom != rec.chrom.shift()) | (rec.bin_start - rec.bin_start.shift() > BIN)).cumsum()
    loci = rec.groupby('locus').agg(chrom=('chrom', 'first'), start=('bin_start', 'min'),
                                    end=('bin_start', 'max'), max_freq=('freq', 'max')).reset_index(drop=True)
    loci['center'] = (loci.start + loci.end + BIN) // 2
    fb = pd.concat([fib[['chrom', 'start']].rename(columns={'start': 'pos'}),
                    fib[['chrom', 'end']].rename(columns={'end': 'pos'})]).sort_values(['chrom', 'pos'])
    fbd = {c: g.pos.values for c, g in fb.groupby('chrom')}

    def nearest(chrom, pos, ref):
        d = np.full(len(pos), np.inf)
        for c in np.unique(chrom):
            r = ref.get(c)
            if r is None or not len(r):
                continue
            sel = chrom == c
            i = np.clip(np.searchsorted(r, pos[sel]), 1, len(r) - 1)
            d[sel] = np.minimum(np.abs(pos[sel] - r[i - 1]), np.abs(pos[sel] - r[i]))
        return d

    loci['dist_fib'] = nearest(loci.chrom.values, loci.center.values, fbd)
    loci['class'] = np.where(loci.dist_fib <= SHARE_DIST, 'shared', 'lcl_specific')
    ld = {c: np.sort(g.center.values) for c, g in loci.groupby('chrom')}
    fb['dist_lcl'] = nearest(fb.chrom.values, fb.pos.values, ld)
    fibspec = fb[fb.dist_lcl > SHARE_DIST].rename(columns={'pos': 'center'})
    fibspec['class'] = 'fib_specific'

    genes = pd.read_csv(GENES, sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'start', 'end', 'name'])
    genes = genes[genes.chrom.isin(nbins)].sort_values(['chrom', 'start'])
    gd = {c: g for c, g in genes.groupby('chrom')}
    gstart = {c: g.start.values for c, g in gd.items()}
    gend_max = {c: np.maximum.accumulate(g.end.values) for c, g in gd.items()}
    allw = pd.concat([loci[['chrom', 'center', 'class', 'max_freq']], fibspec[['chrom', 'center', 'class']]],
                     ignore_index=True)
    allw['genes'] = gene_hits(allw.chrom.values, allw.center.values, gd)
    allw.to_csv(OUT / 'boundary_loci.tsv', sep='\t', index=False)

    rng = np.random.default_rng(1)
    rows = []
    for cls, g in allw.groupby('class'):
        obs = frac_with_gene(g.chrom.values, g.center.values, gstart, gend_max)
        null = []
        for _ in range(N_PERM):
            cc, pp = [], []
            for c, k in g.chrom.value_counts().items():
                cc.append(np.repeat(c, k))
                pp.append(rng.integers(HALF_WIN, nbins[c] * BIN - HALF_WIN, size=k))
            null.append(frac_with_gene(np.concatenate(cc), np.concatenate(pp), gstart, gend_max))
        null = np.array(null)
        rows.append({'class': cls, 'n': len(g), 'frac_with_gene': obs, 'null_mean': null.mean(),
                     'fold': obs / null.mean(), 'p_emp': (np.sum(null >= obs) + 1) / (N_PERM + 1)})
    enr = pd.DataFrame(rows)
    enr.to_csv(OUT / 'boundary_gene_enrichment.tsv', sep='\t', index=False)
    log('\n' + enr.round(4).to_string(index=False))

    gl = (allw.assign(gene=allw.genes.str.split(',')).explode('gene')
          .query('gene != "" and gene == gene').groupby('class').gene.apply(lambda x: sorted(set(x))))
    pd.DataFrame({'class': gl.index, 'n_genes': gl.apply(len), 'genes': gl.apply(','.join)}).to_csv(
        OUT / 'boundary_genes_by_class.tsv', sep='\t', index=False)
    log(gl.apply(len).to_string())
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()


if __name__ == '__main__':
    main()
