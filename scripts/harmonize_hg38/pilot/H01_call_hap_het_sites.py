#!/usr/bin/env python3
"""H01_call_hap_het_sites.py -- call sequence differences between a donor's hap1 and hap2
assemblies, restricted to one hg38 chromosome, as the "het site" evidence for the
read-assignment-confidence filter proposed in md/20260915_haplotype_assignment_modbed_fix.md.

Why this and not statistical/read-based phasing: the two haplotype ASSEMBLIES already
represent the two parental chromosome copies (that's what makes them "hap1"/"hap2" in the
first place) -- so the assembly-vs-assembly difference set is the true, direct haplotype
difference, no BAMs or genotyping needed. Reads whose modbed span doesn't cross any such
difference (see H03) have no information available to have been placed correctly by HPRC2's
score-based assignment (confirmed this session: same forced-tie-split issue as HiFi, see the
proposal doc's §3b resolution) -- they could not have been distinguished even in principle.

Pipeline, both directions (asymmetric alignment means a variant found hap1-as-query isn't
automatically the same call as hap2-as-query, particularly for indels; both are needed to get
proper BED coordinates in EACH haplotype's own coordinate system):
  1. Find the assembly contig(s) covering the target hg38 chromosome via chainmap.py (already
     built + cross-validated against pyliftover this session) on each haplotype's own
     assembly->GRCh38 chain.
  2. Extract those contigs (samtools faidx -- assemblies are already bgzipped+indexed).
  3. minimap2 -cx asm5 --cs, each haplotype as reference in turn.
  4. paftools.js call -> variant BED in the REFERENCE haplotype's coordinates for that run.

Usage: H01_call_hap_het_sites.py <sample> <hg38_chrom>
Writes to <outdir>/<sample>_<chrom>_het_in_hap1.bed and *_het_in_hap2.bed
(outdir = data/het_snps/tmp_<chrom>/<sample>/)
"""
import subprocess
import sys
from pathlib import Path

import chainmap

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
ASSEMBLY_DIR = PROJDIR / 'data' / 'assemblies'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
HET_SNPS_DIR = PROJDIR / 'data' / 'het_snps'

SAMTOOLS = '/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools'
MINIMAP2 = '/u/home/t/terencew/bin/minimap2'
K8 = '/u/home/t/terencew/bin/k8'
PAFTOOLS = '/u/home/t/terencew/bin/paftools.js'


def contigs_for_chrom(sample, hap, chrom):
    chain_path = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'
    blocks = chainmap.parse_chain(str(chain_path))
    return sorted(set(t_name for t_name, t0, t1, q_name, q0, q1, flip in blocks if q_name == chrom))


def extract_contigs(sample, hap, contigs, out_fa):
    fasta = ASSEMBLY_DIR / f'{sample}_hap{hap}.fa.gz'
    names = [f'{sample}#{hap}#{c}' for c in contigs]
    with open(out_fa, 'w') as out:
        subprocess.run([SAMTOOLS, 'faidx', str(fasta)] + names, stdout=out, check=True)


def align_and_call(query_fa, ref_fa, out_paf, out_bed):
    with open(out_paf, 'w') as out:
        subprocess.run([MINIMAP2, '-cx', 'asm5', '--cs', str(ref_fa), str(query_fa)],
                        stdout=out, check=True)
    # paftools.js call emits mixed line types, not plain BED: 'R' lines report the covered
    # reference span (col2=ref_name only, no coordinates in the col2/col3 positions a naive
    # BED parser would expect), 'V' lines are the actual variant calls (col2=ref_name,
    # col3=start, col4=end). Caught before H02 ever ran on this -- a naive
    # `int(line.split('\t')[1])` on an 'R' line raises ValueError on the contig-name string.
    # Keep only 'V' lines, reduced to plain BED3 (contig, start, end) -- H02 only needs
    # interval positions, not the allele/genotype detail also present in the V line.
    call_out = subprocess.run([K8, PAFTOOLS, 'call', str(out_paf)],
                               capture_output=True, text=True, check=True)
    with open(out_bed, 'w') as out:
        for line in call_out.stdout.splitlines():
            f = line.split('\t')
            if f[0] == 'V':
                out.write(f'{f[1]}\t{f[2]}\t{f[3]}\n')


def main():
    sample, chrom = sys.argv[1], sys.argv[2]
    outdir = HET_SNPS_DIR / f'tmp_{chrom}' / sample
    outdir.mkdir(parents=True, exist_ok=True)

    hap1_contigs = contigs_for_chrom(sample, 1, chrom)
    hap2_contigs = contigs_for_chrom(sample, 2, chrom)
    print(f'{sample} {chrom}: hap1 contigs {hap1_contigs}, hap2 contigs {hap2_contigs}',
          flush=True)
    if not hap1_contigs or not hap2_contigs:
        print(f'ERROR: no contigs found for {chrom} on one haplotype', file=sys.stderr)
        sys.exit(1)

    hap1_fa = outdir / 'hap1.fa'
    hap2_fa = outdir / 'hap2.fa'
    extract_contigs(sample, 1, hap1_contigs, hap1_fa)
    extract_contigs(sample, 2, hap2_contigs, hap2_fa)

    # hap2 as reference -> variant calls (paftools.js call's "R" lines, converted to BED by
    # its own logic) are in hap2's OWN coordinates -- what we need for filtering hap2 reads
    align_and_call(hap1_fa, hap2_fa, outdir / 'hap1_vs_hap2.paf', outdir / 'het_in_hap2.bed')
    # and the reverse, for hap1's own coordinates
    align_and_call(hap2_fa, hap1_fa, outdir / 'hap2_vs_hap1.paf', outdir / 'het_in_hap1.bed')

    for hap in (1, 2):
        n = sum(1 for _ in open(outdir / f'het_in_hap{hap}.bed'))
        print(f'{sample} {chrom} hap{hap}: {n} het/difference sites', flush=True)


if __name__ == '__main__':
    main()
