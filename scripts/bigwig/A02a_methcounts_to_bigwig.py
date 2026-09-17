#!/usr/bin/env python3
"""A02a_methcounts_to_bigwig.py -- per-haplotype hg38 methylation bigWig, for visual PMD
sanity-checking in a genome browser (motivation: eyeball whether the large low-methylation
domains the HMM PMD caller found are real, not a calling artifact).

Data source, deliberately: A01a's symmetrized, unfiltered methcounts
(data/pmds/<sample>/<sample>_hap{N}.cpg.methcounts.tsv.gz), NOT P03's production hg38-lifted CpG
matrix. P03 applies H02's k>=1 het-site filter, which drops reads that don't span a
phasing-informative site -- exactly the reads a track meant to show TRUE genome-wide methylation
structure cannot afford to be missing. Using P03's matrix here would show filter-driven gaps as
if they were biology. Same reasoning A01a's own docstring already applied to PMD calling itself;
this script is downstream of the same methcounts, not a fresh reason to filter.

Native assembly-contig coordinates are lifted to hg38 in ONE step at the bedGraph level via
UCSC's own `liftOver` against the existing hap-vs-GRCh38 chain file, rather than lifting each
CpG position individually before aggregation. Chain files are in TARGET=assembly-contig,
QUERY=hg38 orientation (confirmed empirically against a real chain file this session -- liftOver
lifts target->query, i.e. exactly assembly-native -> hg38, using the chain file as-is, no
flipping needed).

Contig-name gotcha (same fix H01/H02 already needed for the modbed/chain mismatch): methcounts
carries the modbed's own SAMPLE#HAP#accession contig naming; the chain file's target name is the
bare accession only. Confirmed directly (first attempt at this script silently produced 0 mapped
lines and 100% "Deleted in new" until this was caught) -- strip the prefix before liftOver.

Usage: A02a_methcounts_to_bigwig.py <sample> <hap> <out.bw>
  Reads data/pmds/<sample>/<sample>_hap<hap>.cpg.methcounts.tsv.gz and
  data/chains/<sample>_hap<hap>_vs_GRCh38.chain.gz.
  Value track is % methylation (0-100), matching common bedMethyl/bigWig conventions.
"""
import gzip
import subprocess
import sys
from pathlib import Path

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
PMD_DIR = PROJDIR / 'data' / 'pmds'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
HG38_CHROMSIZES = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosomal.chromsizes')

KENT = Path('/u/project/cluo/terencew/programs/kentsrc/utils')
LIFTOVER = KENT / 'liftOver'
BEDGRAPH_TO_BW = KENT / 'bedGraphToBigWig'


def methcounts_to_native_bedgraph(methcounts_gz, out_bedgraph):
    """Strip the sample#hap#accession prefix down to the bare accession (chain-file target
    naming), value = % methylation. One row per CpG (methcounts is already symmetrized)."""
    with gzip.open(methcounts_gz, 'rt') as fh, open(out_bedgraph, 'w') as out:
        for line in fh:
            contig, pos, strand, ctx, frac, total = line.rstrip('\n').split('\t')
            bare_contig = contig.rsplit('#', 1)[-1]
            pos = int(pos)
            pct = float(frac) * 100
            out.write(f'{bare_contig}\t{pos}\t{pos + 1}\t{pct:.4f}\n')


def main():
    sample, hap, out_bw = sys.argv[1], sys.argv[2], sys.argv[3]

    methcounts_gz = PMD_DIR / sample / f'{sample}_hap{hap}.cpg.methcounts.tsv.gz'
    chain_gz = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'
    if not methcounts_gz.exists():
        print(f'# SKIP {sample} hap{hap}: no methcounts file', file=sys.stderr)
        sys.exit(0)
    if not chain_gz.exists():
        print(f'# SKIP {sample} hap{hap}: no chain file (known permanent gap for some donors, '
              f'e.g. HG00272)', file=sys.stderr)
        sys.exit(0)

    out_path = Path(out_bw)
    tmpdir = out_path.parent
    native_bg = tmpdir / f'{sample}_hap{hap}.native.bedGraph'
    chain_plain = tmpdir / f'{sample}_hap{hap}.chain'
    hg38_bg = tmpdir / f'{sample}_hap{hap}.hg38.bedGraph'
    unmapped = tmpdir / f'{sample}_hap{hap}.unmapped'
    hg38_sorted_bg = tmpdir / f'{sample}_hap{hap}.hg38.sorted.bedGraph'

    print(f'{sample} hap{hap}: methcounts -> native bedGraph', file=sys.stderr)
    methcounts_to_native_bedgraph(methcounts_gz, native_bg)

    print(f'{sample} hap{hap}: liftOver -> hg38', file=sys.stderr)
    with gzip.open(chain_gz, 'rt') as fh, open(chain_plain, 'w') as out:
        out.write(fh.read())
    subprocess.run([str(LIFTOVER), '-bedPlus=4', str(native_bg), str(chain_plain),
                     str(hg38_bg), str(unmapped)], check=True)

    print(f'{sample} hap{hap}: sort + filter well-formed rows', file=sys.stderr)
    # liftOver's output can include stray malformed rows (caught empirically, not a rare edge
    # case worth ignoring silently) -- require exactly 4 tab-separated fields. Filter+sort as a
    # single shell pipeline (not read-into-Python-then-repipe) -- millions of rows, no reason to
    # hold the whole file in memory twice.
    # Also keep only chroms present in HG38_CHROMSIZES (autosomes -- chrX/chrY/alts would make
    # bedGraphToBigWig abort), and drop repeat starts: two native CpGs can lift onto the same
    # hg38 base, and bedGraphToBigWig rejects overlapping intervals (first one is kept).
    subprocess.run(
        f"awk -F'\\t' 'NR==FNR {{keep[$1]=1; next}} NF==4 && ($1 in keep)' "
        f"{HG38_CHROMSIZES} {hg38_bg} | sort -k1,1 -k2,2n | "
        f"awk -F'\\t' '!($1==c && $2==s) {{print}} {{c=$1; s=$2}}' > {hg38_sorted_bg}",
        shell=True, check=True)

    print(f'{sample} hap{hap}: bedGraphToBigWig', file=sys.stderr)
    subprocess.run([str(BEDGRAPH_TO_BW), str(hg38_sorted_bg), str(HG38_CHROMSIZES),
                     str(out_bw)], check=True)

    for f in (native_bg, chain_plain, hg38_bg, unmapped, hg38_sorted_bg):
        if f.exists():  # allcools env is Python <3.8: no unlink(missing_ok=...)
            f.unlink()

    print(f'# wrote {out_bw}', file=sys.stderr)


if __name__ == '__main__':
    main()
