#!/usr/bin/env python3
"""G01_call_hap_vs_hg38.py -- dipcall-style variant calling of each haplotype assembly against
GRCh38 directly (not hap1-vs-hap2 like H01), restricted to one chromosome for the pilot.

Why this is a SEPARATE product from H01's hap1-vs-hap2 het sites, not a replacement: H01 tells
you where a donor's two haplotypes differ from EACH OTHER, in assembly-native coordinates --
enough to filter modbed reads, but not comparable across donors (each donor's "hap1"/"hap2" is
an arbitrary per-contig label with no cross-donor meaning). This script instead aligns each
haplotype directly to the SHARED GRCh38 reference, giving a phased, hg38-anchored genotype at
every variant site -- the representation needed for cross-donor penetrance/genetic-linkage
analysis (does the SAME hg38 SNP position associate with the SAME methylation direction in
every donor, independent of that donor's own arbitrary hap1/hap2 labeling).

Mirrors asm_lr/scripts/download/D04_run_dipcall.sh's approach (minimap2 asm5 -> paftools.js
call), restricted to one chromosome and reusing this project's already-downloaded assemblies +
reference (GRCh38.autosome.fa, hg38_igvf) rather than a new download. The final merge into one
phased diploid VCF is NOT dipcall-aux.js's vcfpair (tried it first; it requires a GT:AD read-
depth field this assembly-only comparison has no way to populate meaningfully, and threw
"malformatted VCF" on real input as a result -- see merge_diploid()'s docstring) -- phasing is
done directly, which is simple here since it doesn't need to be computed: it comes for free
from the assembly itself (chromosome-scale phased already, no switch-error accumulation over
distance the way statistical/read-based phasing has), the merge just has to encode that
existing phase into one VCF's GT field per shared position.

Usage: G01_call_hap_vs_hg38.py <sample> <hg38_chrom>
Writes to data/vcf_tmp/tmp_<chrom>/<sample>/hap_vs_hg38/:
  hap1.vcf.gz, hap2.vcf.gz    (per-haplotype haploid calls against hg38, sorted+indexed)
  diploid.vcf.gz              (phased 2-haplotype VCF, GT field encodes hap1|hap2 alleles)
"""
import subprocess
import sys
from pathlib import Path

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
ASSEMBLY_DIR = PROJDIR / 'data' / 'assemblies'
VCF_TMP_DIR = PROJDIR / 'data' / 'vcf_tmp'
HERE = Path(__file__).parent

REF_GENOME = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa')

SAMTOOLS = '/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools'
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'
BGZIP = '/u/local/apps/htslib/1.12/gcc-4.8.5/bin/bgzip'
TABIX = '/u/local/apps/htslib/1.12/gcc-4.8.5/bin/tabix'
MINIMAP2 = '/u/home/t/terencew/bin/minimap2'
K8 = '/u/home/t/terencew/bin/k8'
PAFTOOLS = '/u/home/t/terencew/bin/paftools.js'
# dipcall-aux.js (fetched for its 'vcfpair' merge step) turned out not to fit this input shape
# -- see merge_diploid()'s docstring for why the merge is done directly in Python instead.


def extract_ref_chrom(chrom, out_fa):
    if out_fa.exists():
        return
    with open(out_fa, 'w') as out:
        subprocess.run([SAMTOOLS, 'faidx', str(REF_GENOME), chrom], stdout=out, check=True)
    subprocess.run([SAMTOOLS, 'faidx', str(out_fa)], check=True)


def extract_hap_contigs(sample, hap, contigs, out_fa):
    if out_fa.exists():
        return
    fasta = ASSEMBLY_DIR / f'{sample}_hap{hap}.fa.gz'
    names = [f'{sample}#{hap}#{c}' for c in contigs]
    with open(out_fa, 'w') as out:
        subprocess.run([SAMTOOLS, 'faidx', str(fasta)] + names, stdout=out, check=True)


def call_hap_vs_ref(hap_fa, ref_fa, sample_name, outdir):
    paf = outdir / f'{sample_name}.paf'
    sorted_paf = outdir / f'{sample_name}.sorted.paf'
    vcf = outdir / f'{sample_name}.vcf.gz'
    if vcf.exists():
        return vcf

    with open(paf, 'w') as out:
        subprocess.run([MINIMAP2, '-cx', 'asm5', '--cs', str(ref_fa), str(hap_fa)],
                        stdout=out, check=True)
    # paftools.js call requires PAF sorted by target (ref) name, then target start
    subprocess.run(f"sort -k6,6 -k8,8n {paf} > {sorted_paf}", shell=True, check=True)

    raw_vcf = outdir / f'{sample_name}.raw.vcf'
    with open(raw_vcf, 'w') as out:
        subprocess.run([K8, PAFTOOLS, 'call', '-f', str(ref_fa), '-s', sample_name,
                        str(sorted_paf)], stdout=out, check=True)

    subprocess.run(f"{BCFTOOLS} sort {raw_vcf} -Oz -o {vcf}", shell=True, check=True)
    subprocess.run([TABIX, '-p', 'vcf', str(vcf)], check=True)
    return vcf


def _load_hap_variants(vcf_gz):
    """(CHROM, POS) -> (REF, ALT), reading only single-ALT records (paftools.js call's own
    output is always biallelic per record, one ALT per line -- no multi-ALT to normalize)."""
    out = {}
    result = subprocess.run([BCFTOOLS, 'view', '-H', str(vcf_gz)],
                             capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        f = line.split('\t')
        chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
        out[(chrom, pos)] = (ref, alt)
    return out


def _ref_contig_length(chrom, ref_fa):
    fai = Path(str(ref_fa) + '.fai')
    with open(fai) as fh:
        for line in fh:
            f = line.split('\t')
            if f[0] == chrom:
                return int(f[1])
    raise ValueError(f'{chrom} not found in {fai}')


def merge_diploid(hap1_vcf, hap2_vcf, sample, chrom, ref_fa, outdir):
    """Phase hap1's and hap2's independent vs-hg38 calls into one diploid VCF (GT=a|b, hap1
    first) directly in Python, instead of dipcall-aux.js's vcfpair.

    Why not vcfpair: it requires each sample's FORMAT field to be GT:AD (genotype +
    per-allele read depth: `/^(\\.|[0-9]+)\\/(\\.|[0-9]+):(\\S+)/`, checked directly in its
    source, dipcall-aux.js:106) -- a real requirement for dipcall's normal use case (haplotype
    assemblies aligned back to reads for depth support), but meaningless for a pure assembly-
    vs-reference diff, which has no read depth concept at all. paftools.js call's own output
    here is GT only, so vcfpair threw "malformatted VCF" on the very first record fed to it
    (verified: not a formatting slip on this project's end, its parser genuinely requires a
    field this input cannot have). Faking a placeholder AD value would satisfy the regex but
    add fabricated data with no meaning -- writing the phasing logic directly instead, since
    it's simple: each haplotype's call set is already an independent, unambiguous list of
    "this haplotype differs from hg38 here", and phasing them together is just a set union
    keyed by position.
    """
    out_vcf = outdir / 'diploid.vcf.gz'
    if out_vcf.exists():
        return out_vcf

    hap1_vars = _load_hap_variants(hap1_vcf)
    hap2_vars = _load_hap_variants(hap2_vcf)
    all_positions = sorted(set(hap1_vars) | set(hap2_vars))

    plain_vcf = outdir / 'diploid.vcf'
    n_ref_mismatch = 0
    with open(plain_vcf, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write('##source=G01_call_hap_vs_hg38.py (direct phasing, not dipcall-aux.js vcfpair)\n')
        out.write(f'##contig=<ID={chrom},length={_ref_contig_length(chrom, ref_fa)}>\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype, '
                   'hap1|hap2">\n')
        out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n')
        for chrom, pos in all_positions:
            h1 = hap1_vars.get((chrom, pos))
            h2 = hap2_vars.get((chrom, pos))
            if h1 and h2:
                if h1[0] != h2[0]:
                    # same position, different REF called from each haplotype's own alignment
                    # -- an indel-representation mismatch between the two independent
                    # alignments, not a real conflicting genotype. Drop rather than guess;
                    # rare (checked below), and silently picking one REF would misrepresent
                    # the other haplotype's call.
                    n_ref_mismatch += 1
                    continue
                ref = h1[0]
                if h1[1] == h2[1]:
                    alt = h1[1]
                    gt = '1|1'
                else:
                    alt = f'{h1[1]},{h2[1]}'
                    gt = '1|2'
            elif h1:
                ref, alt, gt = h1[0], h1[1], '1|0'
            else:
                ref, alt, gt = h2[0], h2[1], '0|1'
            out.write(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}\n')

    if n_ref_mismatch:
        print(f'{sample}: dropped {n_ref_mismatch} positions with a hap1/hap2 REF-allele '
              f'mismatch (indel-representation disagreement between the two independent '
              f'alignments, not a genuine conflicting genotype)', flush=True)

    subprocess.run(f"{BCFTOOLS} sort {plain_vcf} -Oz -o {out_vcf}", shell=True, check=True)
    subprocess.run([TABIX, '-p', 'vcf', str(out_vcf)], check=True)
    return out_vcf


def main():
    sample, chrom = sys.argv[1], sys.argv[2]
    outdir = VCF_TMP_DIR / f'tmp_{chrom}' / sample / 'hap_vs_hg38'
    outdir.mkdir(parents=True, exist_ok=True)

    # reuse H01's chr15-contig discovery instead of re-deriving it
    sys.path.insert(0, str(HERE))
    import H01_call_hap_het_sites as h01
    hap1_contigs = h01.contigs_for_chrom(sample, 1, chrom)
    hap2_contigs = h01.contigs_for_chrom(sample, 2, chrom)
    print(f'{sample} {chrom}: hap1 contigs {hap1_contigs}, hap2 contigs {hap2_contigs}',
          flush=True)

    ref_fa = outdir / f'{chrom}.fa'
    extract_ref_chrom(chrom, ref_fa)

    hap1_fa = outdir / 'hap1.fa'
    hap2_fa = outdir / 'hap2.fa'
    extract_hap_contigs(sample, 1, hap1_contigs, hap1_fa)
    extract_hap_contigs(sample, 2, hap2_contigs, hap2_fa)

    print('calling hap1 vs hg38...', flush=True)
    hap1_vcf = call_hap_vs_ref(hap1_fa, ref_fa, f'{sample}_hap1', outdir)
    print('calling hap2 vs hg38...', flush=True)
    hap2_vcf = call_hap_vs_ref(hap2_fa, ref_fa, f'{sample}_hap2', outdir)

    print('merging into phased diploid VCF...', flush=True)
    diploid_vcf = merge_diploid(hap1_vcf, hap2_vcf, sample, chrom, ref_fa, outdir)

    n = subprocess.run([BCFTOOLS, 'view', '-H', str(diploid_vcf)],
                        capture_output=True, text=True, check=True).stdout.count('\n')
    print(f'{sample} {chrom}: {n} phased variants in {diploid_vcf}', flush=True)


if __name__ == '__main__':
    main()
