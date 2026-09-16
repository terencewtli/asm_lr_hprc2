#!/usr/bin/env python3
"""QC05_snp_cpg_landscape.py -- validation-scale (10 donors) SNP/CpG landscape QC bundle:
het-SNP counts, CpG counts per haplotype, het-SNP-to-nearest-CpG distance, proportion of reads
covering both a CpG and a het SNP, and CpG-disrupting/creating SNP proportion (using paftools
call's retained .paf files to recover REF/ALT alleles, since H01's het_in_hap{N}.bed only keeps
BED3 positions). chr21 used for the position-level metrics (CpG distance, read-both-cover,
CpG-disrupt/create) for speed; het-SNP counts are genome-wide (cheap, BED text files).

Usage: QC05_snp_cpg_landscape.py
Writes results/qc/data/QC05_*.tsv
"""
import bisect
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
MODBED_DIR = PROJDIR / 'data' / 'modbed'
HET_SNPS_DIR = PROJDIR / 'data' / 'het_snps'
ASSEMBLY_DIR = PROJDIR / 'data' / 'assemblies'
OUTDIR = PROJDIR / 'results' / 'qc' / 'data'
OUTDIR.mkdir(parents=True, exist_ok=True)

SAMTOOLS = '/u/local/apps/samtools/1.15/gcc-4.8.5/bin/samtools'
K8 = '/u/home/t/terencew/bin/k8'
PAFTOOLS = '/u/home/t/terencew/bin/paftools.js'

DONORS = ['HG01884', 'HG01261', 'HG01258', 'HG00438', 'HG00558',
          'HG00290', 'HG00323', 'HG02602', 'HG02698', 'HG01358']
SUPERPOP = {'HG01884': 'AFR', 'HG01261': 'AMR', 'HG01258': 'AMR', 'HG01358': 'AMR',
            'HG00438': 'EAS', 'HG00558': 'EAS', 'HG00290': 'EUR', 'HG00323': 'EUR',
            'HG02602': 'SAS', 'HG02698': 'SAS'}
VALIDATION_CHROM = 'chr21'
ALL_CHROMS = [f'chr{i}' for i in range(1, 23)]


def het_snp_counts():
    rows = []
    for s in DONORS:
        for hap in (1, 2):
            for c in ALL_CHROMS:
                bed = HET_SNPS_DIR / f'tmp_{c}' / s / f'het_in_hap{hap}.bed'
                n = sum(1 for _ in open(bed)) if bed.exists() else None
                rows.append((s, hap, c, n))
    df = pd.DataFrame(rows, columns=['sample', 'hap', 'chrom', 'n_het_snps'])
    df.to_csv(OUTDIR / 'QC05_het_snp_counts.tsv', sep='\t', index=False)
    return df


def het_sites_for_chrom(sample, hap, chrom):
    bed = HET_SNPS_DIR / f'tmp_{chrom}' / sample / f'het_in_hap{hap}.bed'
    sites = defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            sites[f[0]].append(int(f[1]))
    for c in sites:
        sites[c].sort()
    return sites


def nearest_cpg_distance(pos, sorted_cpgs):
    import bisect
    if not sorted_cpgs:
        return None
    i = bisect.bisect_left(sorted_cpgs, pos)
    best = None
    if i < len(sorted_cpgs):
        best = abs(sorted_cpgs[i] - pos)
    if i > 0:
        d = abs(sorted_cpgs[i - 1] - pos)
        if best is None or d < best:
            best = d
    return best


def scan_modbed_region(sample, hap, chrom_contigs, het_sites):
    """Single indexed-fetch (tabix, not a full-file scan) pass over just the contigs relevant
    to VALIDATION_CHROM, computing CpG positions (for the nearest-distance metric) and the
    per-read CpG/het-site coverage counts (for the "covers both" metric) together -- the
    original two-function version did two independent full-file gzip scans per (sample, hap),
    which is what made this painfully slow on a 1-core interactive shell; tabix.fetch() per
    contig is an indexed seek, not a linear scan, and doing both metrics in one pass halves the
    I/O again."""
    modbed_path = MODBED_DIR / f'{sample}_hap{hap}.modbed.gz'
    tbx = pysam.TabixFile(str(modbed_path))
    cpgs = defaultdict(set)
    n_total = n_cpg = n_het = n_both = 0
    for contig in chrom_contigs:
        modbed_contig = f'{sample}#{hap}#{contig}'
        if modbed_contig not in tbx.contigs:
            continue
        sorted_sites = het_sites.get(contig, [])
        for row in tbx.fetch(modbed_contig):
            f = row.split('\t')
            start, end = int(f[1]), int(f[2])
            has_cpg = f[6] not in ('', '.', '-') or f[7] not in ('', '.', '-')
            if has_cpg:
                for field in (f[6], f[7]):
                    if field not in ('', '.', '-'):
                        for off in field.split(','):
                            cpgs[contig].add(start + abs(int(off)))
            lo = bisect.bisect_left(sorted_sites, start)
            hi = bisect.bisect_left(sorted_sites, end)
            has_het = (hi - lo) > 0
            n_total += 1
            n_cpg += has_cpg
            n_het += has_het
            n_both += has_cpg and has_het
    return cpgs, (n_total, n_cpg, n_het, n_both)


def cpg_disrupt_create(sample, contig_seqs):
    """Re-run `k8 paftools.js call` on H01's already-retained .paf files (hap2_vs_hap1.paf ->
    het_in_hap1.bed's source) to recover REF/ALT alleles -- H01's own BED output discards them
    (see H01_call_hap_het_sites.py align_and_call()'s comment). Restricted to SNP-type V lines
    (single-base ref and alt, no indels) on VALIDATION_CHROM's contigs, hap1 side only (ref =
    hap1 assembly, since align_and_call(hap2_fa, hap1_fa, ...) wrote het_in_hap1.bed with hap1
    as paftools' reference).

    contig_seqs: {contig: full_sequence} for hap1, already fetched once by the main loop's
    fetch_contig_seq() calls -- an earlier version of this function spawned one `samtools faidx`
    subprocess PER SNP (tens of thousands per donor on chr21 alone), which was the actual
    bottleneck in this script's first, much-slower run; string-slicing an already-in-memory
    sequence is essentially free by comparison."""
    paf = HET_SNPS_DIR / f'tmp_{VALIDATION_CHROM}' / sample / 'hap2_vs_hap1.paf'
    if not paf.exists():
        return None
    result = subprocess.run([K8, PAFTOOLS, 'call', str(paf)], capture_output=True, text=True)

    snp_records = []
    for line in result.stdout.splitlines():
        f = line.split('\t')
        if f[0] != 'V':
            continue
        contig, start, end, ref, alt = f[1], int(f[2]), int(f[3]), f[6], f[7]
        if len(ref) == 1 and len(alt) == 1 and ref != '-' and alt != '-':
            snp_records.append((contig, start, end, ref, alt))

    n_disrupt = n_create = n_neither = n_skipped = 0
    for contig, start, end, ref, alt in snp_records:
        seq = contig_seqs.get(contig)
        if seq is None or start < 1 or end >= len(seq):
            n_skipped += 1
            continue
        flank = seq[start - 1:end + 1].upper()  # start is 0-based V-line pos; +-1bp
        if len(flank) != 3:
            n_skipped += 1
            continue
        alt_flank = flank[0] + alt.upper() + flank[2]
        ref_has_cpg = 'CG' in flank
        alt_has_cpg = 'CG' in alt_flank
        if ref_has_cpg and not alt_has_cpg:
            n_disrupt += 1
        elif alt_has_cpg and not ref_has_cpg:
            n_create += 1
        else:
            n_neither += 1
    return dict(sample=sample, n_snps_checked=len(snp_records) - n_skipped,
                n_cpg_disrupting=n_disrupt, n_cpg_creating=n_create, n_neither=n_neither)


def fetch_contig_seq(sample, hap, contig):
    """Whole-contig sequence, ONE samtools faidx call (not per-site) -- CpG-architecture
    density metrics (genome-wide and local-around-het-SNP) are sequence features, independent
    of modbed coverage/calls, so this reads directly from the assembly FASTA."""
    fasta = ASSEMBLY_DIR / f'{sample}_hap{hap}.fa.gz'
    region = f'{sample}#{hap}#{contig}'
    out = subprocess.run([SAMTOOLS, 'faidx', str(fasta), region],
                          capture_output=True, text=True, check=True)
    return ''.join(out.stdout.splitlines()[1:]).upper()


def cg_density_per_kb(seq):
    if not seq:
        return None
    return seq.count('CG') / len(seq) * 1000


def local_cg_density(seq, pos, contig_start, window):
    """seq is the FULL contig string (0-based); pos is a genomic (contig-native) coordinate.
    contig_start is always 0 here (whole contig fetched), kept explicit rather than assumed."""
    lo = max(0, pos - window - contig_start)
    hi = min(len(seq), pos + window - contig_start)
    w = seq[lo:hi]
    return cg_density_per_kb(w)


def add_superpop(df, sample_col='sample'):
    df = df.copy()
    df.insert(1, 'superpop', df[sample_col].map(SUPERPOP))
    return df.sort_values(['superpop', sample_col]).reset_index(drop=True)


def write_both_views(df, name, value_cols, agg='mean'):
    """Per-instructions: (1) per-individual view sorted by superpop, (2) superpop-aggregated
    view -- for every metric, not just one aggregate table, so ancestry-level systematic
    differences (a real confound given this project's ancestry-stratified framing) are visible
    alongside per-donor outliers."""
    per_donor = add_superpop(df)
    per_donor.to_csv(OUTDIR / f'QC05_{name}_per_donor.tsv', sep='\t', index=False)
    by_superpop = per_donor.groupby('superpop')[value_cols].agg(agg).reset_index()
    by_superpop.to_csv(OUTDIR / f'QC05_{name}_by_superpop.tsv', sep='\t', index=False)
    return per_donor, by_superpop


def main():
    het_df = het_snp_counts()
    het_per_donor_total = het_df.groupby('sample', as_index=False)['n_het_snps'].sum()
    het_per_donor, het_by_superpop = write_both_views(
        het_per_donor_total, 'het_snp_counts', ['n_het_snps'])
    print('# het SNP counts -- per donor (sorted by superpop):')
    print(het_per_donor)
    print('# het SNP counts -- by superpop (mean total per donor):')
    print(het_by_superpop)

    cpg_rows = []
    dist_rows = []
    both_rows = []
    disrupt_rows = []
    density_rows = []       # per (sample, hap, contig): genome(chr21-scope)-wide CpG/kb
    local_density_rows = [] # per (sample, hap, het site): local CpG/kb at +-100bp / +-500bp

    for s in DONORS:
        hap1_contig_seqs = {}
        for hap in (1, 2):
            het_sites = het_sites_for_chrom(s, hap, VALIDATION_CHROM)
            contigs = set(het_sites.keys())
            cpgs, (n_total, n_cpg, n_het, n_both) = scan_modbed_region(s, hap, contigs, het_sites)
            n_cpg_total = sum(len(v) for v in cpgs.values())
            cpg_rows.append((s, hap, VALIDATION_CHROM, n_cpg_total))

            for contig, sites in het_sites.items():
                sorted_cpgs = sorted(cpgs.get(contig, []))
                for pos in sites:
                    d = nearest_cpg_distance(pos, sorted_cpgs)
                    if d is not None:
                        dist_rows.append((s, hap, d))

                # CpG architecture: sequence-level density, independent of modbed coverage.
                # hap1's sequences are cached for reuse below (cpg_disrupt_create needs hap1
                # as REF, since het_in_hap1.bed's positions are in hap1's own coord frame).
                seq = fetch_contig_seq(s, hap, contig)
                if hap == 1:
                    hap1_contig_seqs[contig] = seq
                density_rows.append((s, hap, contig, cg_density_per_kb(seq), len(seq)))
                for pos in sites:
                    d100 = local_cg_density(seq, pos, 0, 100)
                    d500 = local_cg_density(seq, pos, 0, 500)
                    local_density_rows.append((s, hap, pos, d100, d500))

            both_rows.append((s, hap, n_total, n_cpg, n_het, n_both,
                               n_both / n_total if n_total else None))
            print(f'{s} hap{hap}: done', flush=True)

        dis = cpg_disrupt_create(s, hap1_contig_seqs)
        if dis:
            dis['frac_disrupting'] = dis['n_cpg_disrupting'] / dis['n_snps_checked'] if dis['n_snps_checked'] else None
            dis['frac_creating'] = dis['n_cpg_creating'] / dis['n_snps_checked'] if dis['n_snps_checked'] else None
            disrupt_rows.append(dis)

    cpg_df = pd.DataFrame(cpg_rows, columns=['sample', 'hap', 'chrom', 'n_cpgs'])
    cpg_per_donor_total = cpg_df.groupby('sample', as_index=False)['n_cpgs'].sum()
    cpg_per_donor, cpg_by_superpop = write_both_views(
        cpg_per_donor_total, 'cpg_counts', ['n_cpgs'])

    dist_df = pd.DataFrame(dist_rows, columns=['sample', 'hap', 'het_to_cpg_distance'])
    dist_summary = dist_df.groupby('sample', as_index=False)['het_to_cpg_distance'].median()
    dist_per_donor, dist_by_superpop = write_both_views(
        dist_summary, 'het_cpg_distance_median', ['het_to_cpg_distance'])
    dist_df.to_csv(OUTDIR / 'QC05_het_cpg_distance_raw.tsv', sep='\t', index=False)

    both_df = pd.DataFrame(both_rows, columns=['sample', 'hap', 'n_total_reads',
                                                'n_reads_with_cpg', 'n_reads_with_het',
                                                'n_reads_with_both', 'frac_reads_with_both'])
    both_per_donor_avg = both_df.groupby('sample', as_index=False)['frac_reads_with_both'].mean()
    both_per_donor, both_by_superpop = write_both_views(
        both_per_donor_avg, 'reads_covering_both', ['frac_reads_with_both'])
    both_df.to_csv(OUTDIR / 'QC05_reads_covering_both_raw.tsv', sep='\t', index=False)

    disrupt_df = pd.DataFrame(disrupt_rows)
    disrupt_per_donor, disrupt_by_superpop = write_both_views(
        disrupt_df, 'cpg_disrupt_create', ['frac_disrupting', 'frac_creating'])

    # CpG architecture: genome(chr21-scope)-wide density (README's confound concept -- a
    # sequence feature, computed independent of modbed coverage/calls).
    density_df = pd.DataFrame(density_rows, columns=['sample', 'hap', 'contig',
                                                       'cpg_per_kb', 'contig_len'])
    density_summary = density_df.groupby('sample', as_index=False)['cpg_per_kb'].mean()
    density_per_donor, density_by_superpop = write_both_views(
        density_summary, 'cpg_density_genomewide', ['cpg_per_kb'])
    density_df.to_csv(OUTDIR / 'QC05_cpg_density_genomewide_raw.tsv', sep='\t', index=False)

    # Local CpG density in the window that actually drives k-filter read retention.
    local_df = pd.DataFrame(local_density_rows,
                             columns=['sample', 'hap', 'pos', 'cpg_per_kb_100bp',
                                      'cpg_per_kb_500bp'])
    local_summary = local_df.groupby('sample', as_index=False)[
        ['cpg_per_kb_100bp', 'cpg_per_kb_500bp']].mean()
    local_per_donor, local_by_superpop = write_both_views(
        local_summary, 'cpg_density_local_het_snp', ['cpg_per_kb_100bp', 'cpg_per_kb_500bp'])
    local_df.to_csv(OUTDIR / 'QC05_cpg_density_local_het_snp_raw.tsv', sep='\t', index=False)

    print('\n# CpG counts -- per donor:'); print(cpg_per_donor)
    print('\n# CpG counts -- by superpop:'); print(cpg_by_superpop)
    print('\n# het-to-CpG distance (median bp) -- per donor:'); print(dist_per_donor)
    print('\n# het-to-CpG distance (median bp) -- by superpop:'); print(dist_by_superpop)
    print('\n# overall distance distribution:'); print(dist_df['het_to_cpg_distance'].describe())
    print('\n# reads covering both -- per donor:'); print(both_per_donor)
    print('\n# reads covering both -- by superpop:'); print(both_by_superpop)
    print('\n# CpG disrupt/create -- per donor:'); print(disrupt_per_donor)
    print('\n# CpG disrupt/create -- by superpop:'); print(disrupt_by_superpop)
    print('\n# CpG density genome(chr21-scope)-wide (CpG/kb) -- per donor:'); print(density_per_donor)
    print('\n# CpG density genome-wide -- by superpop:'); print(density_by_superpop)
    print('\n# Local CpG density around het SNPs -- per donor:'); print(local_per_donor)
    print('\n# Local CpG density around het SNPs -- by superpop:'); print(local_by_superpop)


if __name__ == '__main__':
    main()
