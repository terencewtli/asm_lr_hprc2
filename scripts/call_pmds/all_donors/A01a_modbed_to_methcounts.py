#!/usr/bin/env python3
"""A01a_modbed_to_methcounts.py -- replaces the WGBS template's A01a_allc_to_methcounts.sh for
this project's ONT modbed input.

Not a column-rename: modbed is READ-LEVEL (one row per read span, comma-separated per-call
offsets in cols 7/8 -- see docs/data_sources.md Sec.5), not a per-position pileup like allc.
Aggregates every read covering each position into one meth_count/total_count per position,
the same core logic as H02_filtered_locus_matrix.py's parse_locus(), but deliberately WITHOUT
that script's k>=1 het-site filter and WITHOUT liftover to hg38.

Both omissions are deliberate, not shortcuts:
- No het-site filter: P03's production matrix (which DOES apply k>=1) would have been the
  obvious data source to reuse instead of re-deriving this from modbed, but that filter removes
  reads that don't span a phasing-informative site -- exactly the reads a PMD caller cannot
  afford to drop, since it needs comprehensive, unbiased genome-wide coverage. Using P03's
  output here would silently manufacture PMD-like coverage gaps wherever het sites happen to be
  sparse, unrelated to true PMD biology. Confirmed this by reasoning through it before writing
  any code, not caught after the fact.
- No liftover to hg38: PMD calling is a within-haplotype domain-structure question. Working
  directly in each haplotype's own assembly-contig coordinates (matching QC06's approach) avoids
  a liftover step that would only cost precision here for no benefit. Cross-donor/cross-dataset
  PMD comparison, if wanted later, can liftover the resulting PMD intervals (few thousand rows)
  through the existing chain-file infrastructure -- much cheaper than lifting over every CpG
  call up front.

Strand handling and CpG symmetrization (resolved empirically against real data, NOT assumed --
this got corrected mid-implementation, see below): a modbed read's offset sign is fixed by that
READ's own alignment strand (verified: every offset in a `-`-strand read is negative, every
offset in a `+`-strand read is positive -- the sign carries no additional per-call information
beyond col6). Genomic position is `start + abs(offset)`, matching H02's own convention.

First pass at this script assumed (following a surface reading of the WGBS template's A01c,
which appears to call dnmtools pmd on an "unsymmetrized" file) that symmetrization was skippable
-- WRONG, caught by actually reading `dnmtools pmd`'s own usage text: "assumes ... strands are
collapsed so only the positive site appears in the file." The WGBS template's apparent skip was
an illusion: its allc source file is named `*.CGN-Merge.allc.tsv.gz` -- already strand-collapsed
by methylpy upstream of this pipeline entirely, so the template never needed its own merge step.
That precedent does NOT transfer to modbed, which is not pre-merged.

Confirmed directly (not assumed) that `+`-strand calls and `-`-strand calls represent the two
strands of the same CpG dinucleotide at ADJACENT positions, not the same one: in a real 5kb
region, 3803/4255 `+`-strand call positions had a `-`-strand call at position+1, vs. only 164 in
the reverse direction and 17 exact ties -- the textbook signature of `+` = the top-strand C,
`-` = the bottom-strand C reported at the G's coordinate (position+1). Symmetrization therefore
merges a `-`-strand call at position p into position p-1's count, reporting one row per CpG at
the lower (`+`-convention) coordinate, matching dnmtools' documented expectation.

Reference-CpG filter (added 2026-09-16, fixes the "median per-CpG coverage = 1" finding): ~8%
of modbed calls land on positions that are NOT a CG in the haplotype's own assembly (read
errors / mis-assigned reads), each at depth ~1. The first version wrote every called position,
so these junk rows outnumbered real CpGs ~2:1 (NA19338 hap1: 67.8M junk rows at mean depth 1.2
vs 32.2M real CpGs at mean/median depth 30.1/30), dragging the all-row median to 1 and feeding
dnmtools pmd millions of n=1 "CpGs". Only positions where the assembly reads CG at [pos, pos+2)
are now written; dropped rows/calls are logged to stderr per contig.

Usage: A01a_modbed_to_methcounts.py <sample> <hap> <out.methcounts.tsv.gz> [contig ...]
  contig: restrict to specific assembly contigs (e.g. for a chr20-only validation run); omit to
    process every contig in the modbed file.
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pysam

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
MODBED_DIR = PROJDIR / 'data' / 'modbed'
ASSEMBLY_DIR = PROJDIR / 'data' / 'assemblies'


def aggregate_contig(tbx, contig):
    """Symmetrized: a '-'-strand call at position p is folded into position p-1's count (the
    CpG's '+'-convention coordinate), per the empirical strand-adjacency check in this module's
    docstring. Both strands' counts land on the same, single output position."""
    counts = defaultdict(lambda: [0, 0])
    for row in tbx.fetch(contig):
        f = row.split('\t')
        read_start = int(f[1])
        strand = f[5]
        meth_field, unmeth_field = f[6], f[7]
        shift = -1 if strand == '-' else 0
        if meth_field not in ('', '.', '-'):
            for off in meth_field.split(','):
                pos = read_start + abs(int(off)) + shift
                counts[pos][0] += 1
        if unmeth_field not in ('', '.', '-'):
            for off in unmeth_field.split(','):
                pos = read_start + abs(int(off)) + shift
                counts[pos][1] += 1
    return counts


def main():
    sample, hap, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    restrict_contigs = sys.argv[4:] if len(sys.argv) > 4 else None

    modbed_path = MODBED_DIR / f'{sample}_hap{hap}.modbed.gz'
    tbx = pysam.TabixFile(str(modbed_path))
    fasta = pysam.FastaFile(str(ASSEMBLY_DIR / f'{sample}_hap{hap}.fa.gz'))
    contigs = restrict_contigs if restrict_contigs else list(tbx.contigs)

    opener = gzip.open if out_path.endswith('.gz') else open
    n_positions = 0
    n_dropped_rows = n_dropped_calls = 0
    with opener(out_path, 'wt') as out:
        for contig in contigs:
            if contig not in tbx.contigs:
                print(f'# SKIP {sample} hap{hap}: contig {contig} not in modbed', file=sys.stderr)
                continue
            counts = aggregate_contig(tbx, contig)
            seq = fasta.fetch(contig).upper()
            n_kept = drop_rows = drop_calls = 0
            for pos in sorted(counts):
                n_meth, n_unmeth = counts[pos]
                total = n_meth + n_unmeth
                if total == 0:
                    continue
                if seq[pos:pos + 2] != 'CG':
                    drop_rows += 1
                    drop_calls += total
                    continue
                frac = n_meth / total
                out.write(f'{contig}\t{pos}\t+\tCpG\t{frac:.6f}\t{total}\n')
                n_kept += 1
            n_positions += n_kept
            n_dropped_rows += drop_rows
            n_dropped_calls += drop_calls
            print(f'{sample} hap{hap} {contig}: {n_kept} CpG positions kept, '
                  f'{drop_rows} non-CpG positions ({drop_calls} calls) dropped', file=sys.stderr)

    print(f'# wrote {n_positions} methcounts rows to {out_path}; dropped {n_dropped_rows} '
          f'non-CpG positions ({n_dropped_calls} calls)', file=sys.stderr)


if __name__ == '__main__':
    main()
