#!/usr/bin/env python3
"""H02_filtered_locus_matrix.py -- per-CpG methylated/unmethylated counts at a set of hg38
loci, restricted to reads that (a) span >= k het/difference sites between the donor's own
hap1/hap2 assemblies (from H01's output) and (b) don't overlap an HMMFlagger-flagged region.

This is the read-assignment-confidence filter from
md/20260915_haplotype_assignment_modbed_fix.md: HPRC2 assigns every read to a haplotype by
comparing assembly alignment scores and force-splitting exact ties (confirmed this session,
paper Methods text) -- a read with no sequence difference from the other haplotype across its
entire span carries zero evidence for whichever side it landed on. Requiring >=k het sites in
a read's own span is the closest available proxy for "could this read have been distinguished
at all", built without needing the original BAM or read sequence (neither is recoverable from
the modbed alone) -- it does NOT verify the read was assigned to the CORRECT side, only that
assignment was possible in principle (see the proposal doc's §7 for what remains unverified).

Usage: H02_filtered_locus_matrix.py <sample> <region_bed> <k> <out.tsv>
  region_bed: chrom start end name  (hg38 coords, e.g. asm_lr/reference/icr_hg38.bed or a
    Zink-derived bed)
  k: minimum het/difference sites a read's span must contain to be kept (0 = no filter, the
    baseline every other k is compared against)

Requires H01_call_hap_het_sites.py to have already produced
  data/het_snps/tmp_<chrom>/<sample>/het_in_hap{1,2}.bed
for every chrom appearing in region_bed, and data/hmmflagger/<sample>_HMMFlagger.ONT.bed.gz to be
downloaded (see scripts/harmonize_hg38/pilot -- HMMFlagger track fetched directly from the
hprc-epigenome bucket, same pattern as A01a/b/c downloads).

Output: two long TSVs.
  <out.tsv>: per-CpG aggregate -- sample, hap, region, k, hg38_chrom, hg38_pos, n_meth, n_unmeth
  <out.tsv>.reads.tsv: per-READ summary within each locus -- sample, hap, region, k, read_id,
    n_het, n_meth, n_unmeth, frac_meth. This is the one the purity validation (H03) actually
    needs: purity asks whether reads WITHIN one haplotype file agree with each other at a
    locus (a true single haplotype should be internally homogeneous -- reads there either
    mostly-methylated or mostly-unmethylated as a group, not a 50/50 mix), which the per-CpG
    aggregate can't answer since it already collapses across reads.
"""
import bisect
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pysam

import chainmap

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
MODBED_DIR = PROJDIR / 'data' / 'modbed'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
HMM_DIR = PROJDIR / 'data' / 'hmmflagger'
HET_SNPS_DIR = PROJDIR / 'data' / 'het_snps'


def load_het_sites(sample, chrom, hap):
    """contig (bare accession, no SAMPLE#hap# prefix) -> sorted list of het/difference
    positions, for bisect-based range counting.

    H01's het_in_hap{N}.bed uses paftools.js call's own contig naming, which is the FULL
    FASTA header used for that alignment ("SAMPLE#hap#accession", since that's the assembly's
    own header format) -- NOT the bare accession chainmap.q_interval_to_t() returns (chain
    files use bare accessions as t_name). Caught before this ever ran for real: a first version
    of this function kept the full prefixed name as the dict key, which the parse_locus lookup
    below (keyed by chainmap's bare accession) would never match -- every read's n_het would
    have silently come out 0, dropping every read at any k>=1 with no error. Strip the prefix
    here so both sides agree, same fix pattern already used for the modbed/chain mismatch in
    P02_liftover_region_to_hg38.py.
    """
    path = HET_SNPS_DIR / f'tmp_{chrom}' / sample / f'het_in_hap{hap}.bed'
    sites = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split('\t')
            contig, start = f[0].rsplit('#', 1)[-1], int(f[1])
            sites[contig].append(start)
    for contig in sites:
        sites[contig].sort()
    return sites


def load_hmmflagger(sample):
    path = HMM_DIR / f'{sample}_HMMFlagger.ONT.bed.gz'
    intervals = defaultdict(list)
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            f = line.split('\t')
            contig, start, end = f[0], int(f[1]), int(f[2])
            intervals[contig].append((start, end))
    for contig in intervals:
        intervals[contig].sort()
    return intervals


def overlaps_flagged(intervals_for_contig, start, end):
    # linear scan is fine -- HMMFlagger tracks are sparse (hundreds of intervals per contig,
    # not millions), no need for an interval tree here
    for s, e in intervals_for_contig:
        if s < end and start < e:
            return True
        if s >= end:
            break
    return False


def het_count_in_span(sorted_sites, start, end):
    lo = bisect.bisect_left(sorted_sites, start)
    hi = bisect.bisect_left(sorted_sites, end)
    return hi - lo


_CHAINMAP_CACHE = {}
_TABIX_CACHE = {}


def get_chainmap(sample, hap):
    """ChainMap() parses the whole chain file (tens of thousands of lines) on construction --
    caching by (sample, hap) turns an accidental O(n_calls) re-parse into O(n_sample_haps).
    Caught by H03's k-sweep taking >25 min and still climbing for a single donor/chromosome
    (98 loci x 2 haps x 5 k values = 980 parse_locus calls, each rebuilding the parser from
    scratch) -- killed and fixed rather than just waiting it out, since every future
    (more loci, more donors, whole-chromosome) run would hit the same wall harder."""
    key = (sample, hap)
    if key not in _CHAINMAP_CACHE:
        chain_path = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'
        _CHAINMAP_CACHE[key] = chainmap.ChainMap(str(chain_path))
    return _CHAINMAP_CACHE[key]


def get_tabix(sample, hap):
    key = (sample, hap)
    if key not in _TABIX_CACHE:
        modbed_path = MODBED_DIR / f'{sample}_hap{hap}.modbed.gz'
        _TABIX_CACHE[key] = pysam.TabixFile(str(modbed_path))
    return _TABIX_CACHE[key]


def parse_locus(sample, hap, region_name, chrom, start, end, het_sites, hmm_intervals, k):
    cm = get_chainmap(sample, hap)
    hits = cm.q_interval_to_t(chrom, start, end)
    if not hits:
        return [], [], 0, 0

    rows = []
    read_rows = []
    n_total = n_kept = 0
    tbx = get_tabix(sample, hap)

    for asm_contig, asm_start, asm_end in hits:
        modbed_contig = f'{sample}#{hap}#{asm_contig}'
        sorted_sites = het_sites.get(asm_contig, [])
        hmm_for_contig = hmm_intervals.get(modbed_contig, [])
        counts = defaultdict(lambda: [0, 0])

        for row in tbx.fetch(modbed_contig, asm_start, asm_end):
            f = row.split('\t')
            read_id = f[3]
            read_start, read_end = int(f[1]), int(f[2])
            n_total += 1

            n_het = het_count_in_span(sorted_sites, read_start, read_end)
            if n_het < k:
                continue
            if overlaps_flagged(hmm_for_contig, read_start, read_end):
                continue
            n_kept += 1

            read_n_meth = read_n_unmeth = 0
            meth_field, unmeth_field = f[6], f[7]
            if meth_field not in ('', '.', '-'):
                for off in meth_field.split(','):
                    pos = read_start + abs(int(off))
                    if asm_start <= pos < asm_end:
                        counts[pos][0] += 1
                        read_n_meth += 1
            if unmeth_field not in ('', '.', '-'):
                for off in unmeth_field.split(','):
                    pos = read_start + abs(int(off))
                    if asm_start <= pos < asm_end:
                        counts[pos][1] += 1
                        read_n_unmeth += 1

            if read_n_meth + read_n_unmeth > 0:
                frac = read_n_meth / (read_n_meth + read_n_unmeth)
                read_rows.append((sample, hap, region_name, k, read_id, n_het,
                                   read_n_meth, read_n_unmeth, frac))

        for pos in sorted(counts):
            mapped = cm.t_to_q(asm_contig, pos)
            if mapped is None:
                continue
            hg38_chrom, hg38_pos = mapped
            n_meth, n_unmeth = counts[pos]
            rows.append((sample, hap, region_name, k, hg38_chrom, hg38_pos, n_meth, n_unmeth))

    return rows, read_rows, n_total, n_kept


def main():
    sample, region_bed, k, out_path = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]

    het_sites_cache = {}
    hmm_intervals = load_hmmflagger(sample)

    all_rows = []
    all_read_rows = []
    with open(region_bed) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            chrom, start, end, name = f[0], int(f[1]), int(f[2]), f[3]

            for hap in (1, 2):
                cache_key = (chrom, hap)
                if cache_key not in het_sites_cache:
                    het_sites_cache[cache_key] = load_het_sites(sample, chrom, hap)
                rows, read_rows, n_total, n_kept = parse_locus(
                    sample, hap, name, chrom, start, end,
                    het_sites_cache[cache_key], hmm_intervals, k)
                all_rows.extend(rows)
                all_read_rows.extend(read_rows)
                frac = n_kept / n_total if n_total else float('nan')
                print(f'{sample} hap{hap} {name} k={k}: {n_kept}/{n_total} reads kept '
                      f'({frac:.2%}), {len(rows)} CpGs', file=sys.stderr)

    with open(out_path, 'w') as out:
        out.write('sample\thap\tregion\tk\thg38_chrom\thg38_pos\tn_meth\tn_unmeth\n')
        for row in all_rows:
            out.write('\t'.join(str(x) for x in row) + '\n')
    print(f'# wrote {len(all_rows)} rows to {out_path}', file=sys.stderr)

    reads_path = out_path + '.reads.tsv'
    with open(reads_path, 'w') as out:
        out.write('sample\thap\tregion\tk\tread_id\tn_het\tn_meth\tn_unmeth\tfrac_meth\n')
        for row in all_read_rows:
            out.write('\t'.join(str(x) for x in row) + '\n')
    print(f'# wrote {len(all_read_rows)} rows to {reads_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
