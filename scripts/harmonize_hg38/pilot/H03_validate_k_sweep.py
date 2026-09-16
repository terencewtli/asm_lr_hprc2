#!/usr/bin/env python3
"""H03_validate_k_sweep.py -- the actual validation step from
md/20260915_haplotype_assignment_modbed_fix.md §4: sweep k over the Zink et al. 2018
parent-of-origin DMRs (98 on chr15, asm_lr/data/zink_2018_supp5_pofo_dmrs.csv, confirmed hg38
per asm_lr/md/pipeline.md) and check whether within-haplotype read purity rises and plateaus
as k increases -- the direct test of whether the het-site filter is doing what it's supposed
to (removing reads that couldn't have been correctly assigned) rather than just removing reads
for no reason.

Purity definition (per region x haplotype x k): reads within ONE haplotype file, at a locus
that is genuinely imprinted, should be internally homogeneous (either essentially all-
methylated or all-unmethylated as a group) if the haplotype split is correct -- a mix of both
within the "same" haplotype at a truly bimodal locus is the signature of contamination from
force-assigned, uninformative reads. purity = max(n_high, n_low) / n_reads, where a read is
"high" if its own frac_meth >= 0.5. Averaged across all (region, hap) pairs with >= MIN_READS
reads at that k (too few reads makes purity a coin flip, not a real measurement).

Usage: H03_validate_k_sweep.py <sample> <region_bed> <out_prefix>
Writes <out_prefix>_summary.tsv (one row per k: mean purity, n loci included, total reads
retained) and <out_prefix>_per_locus.tsv (per region x hap x k, for inspection).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from H02_filtered_locus_matrix import load_het_sites, load_hmmflagger, parse_locus

K_SWEEP = [0, 1, 2, 3, 5, 10, 20, 30]
MIN_READS = 3


def purity(read_rows):
    if len(read_rows) < MIN_READS:
        return None
    n_high = sum(1 for r in read_rows if r[-1] >= 0.5)
    n_low = len(read_rows) - n_high
    return max(n_high, n_low) / len(read_rows)


def main():
    sample, region_bed, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

    regions = []
    with open(region_bed) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            regions.append((f[0], int(f[1]), int(f[2]), f[3]))

    hmm_intervals = load_hmmflagger(sample)
    het_sites_cache = {}

    per_locus_rows = []
    summary_rows = []

    for k in K_SWEEP:
        purities = []
        total_kept = total_seen = 0
        for chrom, start, end, name in regions:
            for hap in (1, 2):
                cache_key = (chrom, hap)
                if cache_key not in het_sites_cache:
                    het_sites_cache[cache_key] = load_het_sites(sample, chrom, hap)
                _, read_rows, n_total, n_kept = parse_locus(
                    sample, hap, name, chrom, start, end,
                    het_sites_cache[cache_key], hmm_intervals, k)
                total_kept += n_kept
                total_seen += n_total
                p = purity(read_rows)
                per_locus_rows.append((sample, name, hap, k, n_total, n_kept,
                                        len(read_rows), p if p is not None else ''))
                if p is not None:
                    purities.append(p)

        mean_purity = sum(purities) / len(purities) if purities else float('nan')
        retention = total_kept / total_seen if total_seen else float('nan')
        summary_rows.append((k, mean_purity, len(purities), total_kept, total_seen, retention))
        print(f'k={k}: mean_purity={mean_purity:.4f} over {len(purities)} (locus,hap) pairs '
              f'with >={MIN_READS} reads; retention {total_kept}/{total_seen} ({retention:.2%})',
              flush=True)

    with open(f'{out_prefix}_summary.tsv', 'w') as out:
        out.write('k\tmean_purity\tn_loci_with_enough_reads\ttotal_reads_kept\t'
                   'total_reads_seen\tretention_frac\n')
        for row in summary_rows:
            out.write('\t'.join(str(x) for x in row) + '\n')

    with open(f'{out_prefix}_per_locus.tsv', 'w') as out:
        out.write('sample\tregion\thap\tk\tn_reads_total\tn_reads_kept\t'
                   'n_reads_with_calls\tpurity\n')
        for row in per_locus_rows:
            out.write('\t'.join(str(x) for x in row) + '\n')

    print(f'wrote {out_prefix}_summary.tsv and {out_prefix}_per_locus.tsv', flush=True)


if __name__ == '__main__':
    main()
