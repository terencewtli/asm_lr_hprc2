#!/usr/bin/env python3
"""P01_parse_modbed_region.py -- parse a modbed region into per-CpG methylated/unmethylated
counts, in ASSEMBLY (haplotype) coordinates.

Decoding per docs/data_sources.md §5 (confirmed against modbedtools source, not just the
published spec, which the earlier session found was wrong about the offset sign):
    col1 chrom (assembly contig, "SAMPLE#hap#accession")
    col2 start (0-based)          col3 end
    col4 read_id  col5 score  col6 strand
    col7 comma-separated METHYLATED-base offsets, relative to col2
    col8 comma-separated UNMETHYLATED-base offsets, relative to col2
    genomic_position = col2 + abs(offset); sign encodes strand/stored-base of the call, not
    used here (haplotype-assembly coordinate is already strand-fixed by contig orientation --
    flagged in data_sources.md as "worth confirming with one real worked example before
    assuming", which P02's round-trip + bigwig cross-check does).

Usage: P01_parse_modbed_region.py <modbed.gz> <contig> <start> <end>
Prints a TSV to stdout: contig  pos  n_meth  n_unmeth   (pos is 0-based, assembly coordinate)
"""
import sys
from collections import defaultdict

import pysam


def parse_region(modbed_path, contig, start, end):
    """tabix .fetch(contig, start, end) returns any READ whose span overlaps the window -- ONT
    reads are tens to hundreds of kb, so a fetched read's own calls extend far beyond the
    queried window in both directions. Caught empirically (a first pass without the `start <=
    pos < end` filter below produced positions >100kb outside a 3.4kb query, one read at a
    time): every decoded position must be clipped to the requested region explicitly, fetching
    reads is not the same as fetching positions."""
    tbx = pysam.TabixFile(modbed_path)
    counts = defaultdict(lambda: [0, 0])  # pos -> [n_meth, n_unmeth]
    n_reads = 0
    for row in tbx.fetch(contig, start, end):
        f = row.split('\t')
        read_start = int(f[1])
        meth_field, unmeth_field = f[6], f[7]
        if meth_field not in ('', '.', '-'):
            for off in meth_field.split(','):
                pos = read_start + abs(int(off))
                if start <= pos < end:
                    counts[pos][0] += 1
        if unmeth_field not in ('', '.', '-'):
            for off in unmeth_field.split(','):
                pos = read_start + abs(int(off))
                if start <= pos < end:
                    counts[pos][1] += 1
        n_reads += 1
    return counts, n_reads


def main():
    modbed_path, contig, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    counts, n_reads = parse_region(modbed_path, contig, start, end)
    print(f'# {n_reads} reads, {len(counts)} distinct positions with a call', file=sys.stderr)
    for pos in sorted(counts):
        n_meth, n_unmeth = counts[pos]
        print(f'{contig}\t{pos}\t{n_meth}\t{n_unmeth}')


if __name__ == '__main__':
    main()
