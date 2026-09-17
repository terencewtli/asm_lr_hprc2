#!/usr/bin/env python3
"""C02a_call_asm.py -- read-level ASM calls for one (sample, hg38 chromosome), plus the cleaned
het-filtered per-CpG counts. Replaces P03_build_donor_chrom_matrix.py.

Why P03 is superseded (checked 2026-09-17): P03 reused H02.parse_locus, which
  (1) keeps calls at any position (no assembly-CpG check; the junk-row bug fixed in A01a),
  (2) doesn't merge strands (a '-' read's call sits at the CpG's G, one base right),
  (3) lifts with ChainMap.t_to_q, which has no strand-flip correction, and
  (4) throws away read identity, which the read-level test needs.
It keeps H02's read filters unchanged: >= K het/difference sites inside the read span (H01) and
no overlap with an HMMFlagger-flagged assembly interval.

Per haplotype:
  - modbed reads are fetched for every assembly contig with chain blocks on <chrom>
  - reads failing the het/HMMFlagger filters are dropped
  - each call: pos = start + |offset|, minus 1 for '-' strand reads; kept only if the assembly
    reads CG there; mapped to hg38 with meth_bins/B01a_hap_bins.map_points; kept if on <chrom>
Then, for every region in results/asm/regions/regions_hg38.tsv.gz on <chrom>:
  - each read's fraction methylated over the region (reads with >= MIN_CALLS calls there)
  - ont_asm_caller.test_region_reads(hap1 fractions, hap2 fractions): Welch t / Mann-Whitney
    (conservative of the two), or exact separation; needs >= MIN_READS reads per hap

Outputs (results/asm/):
  calls/<sample>_<chrom>.asm.tsv.gz   one row per tested region
  cpg/<sample>_<chrom>.cpg.tsv.gz     hap, pos, n_meth, n_unmeth (het-filtered, hg38)
  calls/<sample>_<chrom>.log.tsv      read/call filter counts per hap

Usage: C02a_call_asm.py <sample> <chrom> [k]
"""
import sys
from pathlib import Path

import math

import numpy as np
import pandas as pd
import pysam
from scipy.special import comb

if not hasattr(math, 'comb'):  # allcools env is Python 3.7; ont_asm_caller uses math.comb (3.8+)
    math.comb = lambda n, k: comb(n, k, exact=True)

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
sys.path.insert(0, str(PROJDIR / 'scripts' / 'meth_bins'))
sys.path.insert(0, str(PROJDIR / 'scripts' / 'harmonize_hg38' / 'pilot'))
sys.path.insert(0, str(PROJDIR / 'github' / 'ont_asm_caller'))
from B01a_hap_bins import load_blocks, map_points  # noqa: E402
from H02_filtered_locus_matrix import load_het_sites, load_hmmflagger  # noqa: E402
from ont_asm_caller import test_region_reads  # noqa: E402

MODBED_DIR = PROJDIR / 'data' / 'modbed'
ASM_DIR = PROJDIR / 'data' / 'assemblies'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
REGIONS = PROJDIR / 'results' / 'asm' / 'regions' / 'regions_hg38.tsv.gz'
OUT = PROJDIR / 'results' / 'asm'
DEFAULT_K = 1
MIN_CALLS = 3
MIN_READS = 3
RANGE_GAP = 1_000_000


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def offsets(field: str) -> np.ndarray:
    if field in ('', '.', '-'):
        return np.empty(0, dtype=np.int64)
    return np.abs(np.array(field.split(','), dtype=np.int64))


def fetch_ranges(tbx, full: str, b: pd.DataFrame, chrom: str):
    """Only the assembly ranges whose chain blocks land on <chrom> (merged within RANGE_GAP),
    so a contig with a small piece on <chrom> isn't scanned end to end."""
    b = b[b.q == chrom].sort_values('t0')
    s, e = b.t0.values, np.maximum.accumulate(b.t1.values)
    brk = np.r_[True, s[1:] > e[:-1] + RANGE_GAP]
    for lo, hi in zip(s[brk], np.r_[e[np.flatnonzero(brk)[1:] - 1], e[-1]]):
        yield from tbx.fetch(full, int(lo), int(hi))


def hap_calls(sample: str, hap: int, chrom: str, k: int, hmm: dict) -> tuple:
    blocks = load_blocks(CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz')
    contigs = [t for t, b in blocks.items() if (b.q == chrom).any()]
    het = load_het_sites(sample, chrom, hap)
    tbx = pysam.TabixFile(str(MODBED_DIR / f'{sample}_hap{hap}.modbed.gz'))
    fa = pysam.FastaFile(str(ASM_DIR / f'{sample}_hap{hap}.fa.gz'))
    stats = dict(sample=sample, hap=hap, chrom=chrom, contigs=len(contigs), reads=0,
                 reads_het=0, reads_kept=0, calls=0, calls_cpg=0, calls_mapped=0)
    rid_l, pos_l, call_l = [], [], []
    rid = 0
    for bare in contigs:
        full = f'{sample}#{hap}#{bare}'
        if full not in tbx.contigs:
            continue
        seq = np.frombuffer(fa.fetch(full).upper().encode(), dtype=np.uint8)
        is_cg = np.zeros(len(seq), dtype=bool)
        is_cg[:-1] = (seq[:-1] == ord('C')) & (seq[1:] == ord('G'))
        sites = np.asarray(het.get(bare, []), dtype=np.int64)
        flagged = np.asarray(hmm.get(full, []), dtype=np.int64).reshape(-1, 2)
        fl_s, fl_e = flagged[:, 0], np.maximum.accumulate(flagged[:, 1]) if len(flagged) else flagged[:, 1]
        c_pos, c_call, c_rid = [], [], []
        seen = set()
        for row in fetch_ranges(tbx, full, blocks[bare], chrom):
            f = row.split('\t')
            if f[3] in seen:
                continue
            seen.add(f[3])
            s, e = int(f[1]), int(f[2])
            stats['reads'] += 1
            if np.searchsorted(sites, e) - np.searchsorted(sites, s) < k:
                continue
            stats['reads_het'] += 1
            if len(fl_s):
                i = np.searchsorted(fl_s, e) - 1   # intervals starting before read end
                if i >= 0 and fl_e[i] > s:
                    continue
            stats['reads_kept'] += 1
            shift = -1 if f[5] == '-' else 0
            om, ou = offsets(f[6]), offsets(f[7])
            p = np.concatenate([om, ou]) + s + shift
            c = np.concatenate([np.ones(len(om), np.int8), np.zeros(len(ou), np.int8)])
            c_pos.append(p)
            c_call.append(c)
            c_rid.append(np.full(len(p), rid, dtype=np.int32))
            rid += 1
        if not c_pos:
            continue
        p = np.concatenate(c_pos)
        c = np.concatenate(c_call)
        r = np.concatenate(c_rid)
        stats['calls'] += len(p)
        ok = (p >= 0) & (p < len(seq) - 1)
        ok[ok] = is_cg[p[ok]]
        p, c, r = p[ok], c[ok], r[ok]
        stats['calls_cpg'] += len(p)
        mok, qn, qp = map_points(blocks[bare], p)
        keep = mok & (qn == chrom)
        stats['calls_mapped'] += int(keep.sum())
        pos_l.append(qp[keep])
        call_l.append(c[keep])
        rid_l.append(r[keep])
    if not pos_l:
        return pd.DataFrame(columns=['rid', 'pos', 'call']), stats
    df = pd.DataFrame({'rid': np.concatenate(rid_l), 'pos': np.concatenate(pos_l),
                       'call': np.concatenate(call_l)})
    return df.drop_duplicates(['rid', 'pos']), stats


def read_fracs(df: pd.DataFrame, reg: pd.DataFrame) -> pd.DataFrame:
    """(region index, read) -> fraction methylated, for reads with >= MIN_CALLS calls."""
    starts, ends = reg.start.values, reg.end.values
    i = np.searchsorted(starts, df.pos.values, side='right') - 1
    inside = (i >= 0) & (df.pos.values < ends[np.clip(i, 0, None)])
    d = pd.DataFrame({'reg': i[inside], 'rid': df.rid.values[inside], 'call': df.call.values[inside]})
    g = d.groupby(['reg', 'rid']).call.agg(['sum', 'size']).reset_index()
    g = g[g['size'] >= MIN_CALLS]
    g['frac'] = g['sum'] / g['size']
    return g


def main() -> None:
    sample, chrom = sys.argv[1], sys.argv[2]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_K
    tag = f'{sample}_{chrom}'
    for d in ('calls', 'cpg'):
        (OUT / d).mkdir(parents=True, exist_ok=True)

    hmm = load_hmmflagger(sample)
    calls, stats = {}, []
    for hap in (1, 2):
        try:
            calls[hap], st = hap_calls(sample, hap, chrom, k, hmm)
        except FileNotFoundError as err:
            log(f'# SKIP {tag} hap{hap}: {err}')
            continue
        stats.append(st)
        log(f'{tag} hap{hap}: ' + ' '.join(f'{a}={b}' for a, b in st.items() if a not in ('sample', 'hap', 'chrom')))
    pd.DataFrame(stats).to_csv(OUT / 'calls' / f'{tag}.log.tsv', sep='\t', index=False)

    cpg = []
    for hap, df in calls.items():
        g = df.groupby('pos').call.agg(['sum', 'size']).reset_index()
        cpg.append(pd.DataFrame({'hap': hap, 'pos': g.pos, 'n_meth': g['sum'],
                                 'n_unmeth': g['size'] - g['sum']}))
    if cpg:
        pd.concat(cpg).to_csv(OUT / 'cpg' / f'{tag}.cpg.tsv.gz', sep='\t', index=False)

    cols = ['region_id', 'set', 'chrom', 'start', 'end', 'n_cpg_ref', 'n_reads1', 'n_reads2',
            'mu1', 'mu2', 'delta', 'stat', 'pval', 'method']
    rows = []
    if len(calls) == 2:
        reg = pd.read_csv(REGIONS, sep='\t')
        reg = reg[reg.chrom == chrom]
        for rset, r in reg.groupby('set'):
            r = r.sort_values('start').reset_index(drop=True)
            fr = {h: read_fracs(calls[h], r) for h in (1, 2)}
            by = {h: {i: g.frac.values for i, g in fr[h].groupby('reg')} for h in (1, 2)}
            for i in sorted(set(by[1]) & set(by[2])):
                f1, f2 = by[1][i], by[2][i]
                if len(f1) < MIN_READS or len(f2) < MIN_READS:
                    continue
                res = test_region_reads(chrom, int(r.start[i]), int(r.end[i]), f1[:, None],
                                        f2[:, None], min_calls=1, min_reads=MIN_READS)
                if res is None:
                    continue
                rows.append((r.region_id[i], rset, chrom, res.start, res.end, int(r.n_cpg_ref[i]),
                             res.n_reads1, res.n_reads2, res.mu1_hat, res.mu2_hat, res.delta,
                             res.stat, res.pval, res.method))
            log(f'{tag} {rset}: {len(r):,} regions, tested so far {len(rows):,}')
    out = pd.DataFrame(rows, columns=cols)
    out.insert(0, 'sample', sample)
    out.to_csv(OUT / 'calls' / f'{tag}.asm.tsv.gz', sep='\t', index=False, float_format='%.6g')
    log(f'# wrote {len(out):,} tested regions for {tag}')


if __name__ == '__main__':
    main()
