#!/usr/bin/env python3
"""P02_liftover_region_to_hg38.py -- project P01's per-CpG assembly-coordinate counts to hg38
via chainmap.py, keeping BOTH coordinate systems (per the "don't discard the original haplotype
coordinate" note in md/20260915.chatgpt.md).

Usage: P02_liftover_region_to_hg38.py <p01_output.tsv> <chain.gz>
Prints a TSV to stdout: asm_contig  asm_pos  hg38_chrom  hg38_pos  n_meth  n_unmeth
Positions that fail to lift (fall in a chain gap -- indel/rearrangement) are dropped, counted,
and reported on stderr rather than silently omitted.
"""
import sys

import chainmap


def main():
    p01_path, chain_path = sys.argv[1], sys.argv[2]
    cm = chainmap.ChainMap(chain_path)

    n_total = n_lifted = 0
    rows = []
    with open(p01_path) as fh:
        for line in fh:
            contig, pos, n_meth, n_unmeth = line.rstrip('\n').split('\t')
            pos = int(pos)
            n_total += 1
            # modbed contigs are "SAMPLE#hap#accession" (e.g. "HG00126#1#CM090115.1"); chain
            # files' t_name is the bare accession ("CM090115.1") -- caught empirically (a first
            # pass silently lifted 0/399 positions here) rather than assumed from the docstring
            # note in chainmap.py.
            bare_contig = contig.rsplit('#', 1)[-1]
            mapped = cm.t_to_q(bare_contig, pos)
            if mapped is None:
                continue
            n_lifted += 1
            hg38_chrom, hg38_pos = mapped
            rows.append((contig, pos, hg38_chrom, hg38_pos, n_meth, n_unmeth))

    print(f'# {n_lifted}/{n_total} positions lifted '
          f'({n_total - n_lifted} fell in a chain gap)', file=sys.stderr)

    # cross-check a sample against pyliftover before trusting the bulk output
    if rows:
        sample_contig = rows[0][0].rsplit('#', 1)[-1]
        sample_positions = [r[1] for r in rows if r[0].rsplit('#', 1)[-1] == sample_contig]
        checked, agree, bad = chainmap.validate_against_pyliftover(
            chain_path, cm, sample_contig, sample_positions, n=30)
        print(f'# pyliftover cross-check: {agree}/{checked} agree', file=sys.stderr)
        if bad:
            print(f'# DISAGREEMENTS: {bad[:5]}', file=sys.stderr)

    rows.sort(key=lambda r: (r[2], r[3]))
    for contig, pos, hg38_chrom, hg38_pos, n_meth, n_unmeth in rows:
        print(f'{contig}\t{pos}\t{hg38_chrom}\t{hg38_pos}\t{n_meth}\t{n_unmeth}')


if __name__ == '__main__':
    main()
