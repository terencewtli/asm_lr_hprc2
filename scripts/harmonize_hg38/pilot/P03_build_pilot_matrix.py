#!/usr/bin/env python3
"""P03_build_pilot_matrix.py -- run H02's het-filtered, liftover'd locus parser for ONE sample
(both haplotypes) across a set of hg38 regions, producing that donor's slice of the eventual
donor x CpG table.

Rewritten (2026-09-15) to call H02_filtered_locus_matrix directly instead of orchestrating raw
P01+P02 -- H02 is a strict superset of what P01+P02 did (parses modbed, lifts over via
chainmap, AND applies the het-site + HMMFlagger filter), so duplicating that logic here instead
of reusing it would just be two copies to keep in sync. The ORIGINAL unfiltered 229-donor run
from this script's prior version is deliberately left in place at results/pilot/per_sample/ --
not overwritten, not deleted -- as a baseline to compare the filtered output against, not
superseded silently.

k chosen (see md/20260915_haplotype_assignment_modbed_fix.md and this session's follow-up
validation): default k=1, not the doc's originally-proposed k>=2-3. The k=0..30 sweep against
98 chr15 Zink DMRs (H03) showed purity essentially flat across that whole range (0.7545 at k=0
to 0.7568 at k=30) -- filtering harder than k=1 buys no measured benefit, so there's no
empirical case for paying the extra retention cost. k=1 is still worth keeping as a floor: a
read crossing zero het sites has literally zero assignment evidence regardless of what it does
to any aggregate metric, and it's nearly free (drops <0.5% of reads in the pilot). Separately,
GMM classification (H04) showed only ~33% of loci are even genuinely bimodal -- the other ~67%
have no real two-state structure for a phasing problem to manifest in at all, which is the
right context for why the k-sweep found nothing: most of this pilot's loci aren't the kind of
locus this filter could show an effect at in the first place.

Usage: P03_build_pilot_matrix.py <sample_id> <region_bed> <out.tsv> [k]
  region_bed: chrom start end name  (e.g. asm_lr/reference/icr_hg38.bed)
  k: minimum het/difference sites a read's span must contain to be kept (default 1)

Output: long-format TSV -- sample, hap, region_name, hg38_chrom, hg38_pos, n_meth, n_unmeth,
frac_meth -- same shape as the original unfiltered version, so downstream code that reads
either is unaffected by which one it's pointed at.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from H02_filtered_locus_matrix import load_hmmflagger, load_het_sites, parse_locus

DEFAULT_K = 1


def main():
    sample, region_bed, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    k = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_K

    hmm_intervals = load_hmmflagger(sample)
    het_sites_cache = {}

    all_rows = []
    with open(region_bed) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            chrom, start, end, name = f[0], int(f[1]), int(f[2]), f[3]

            for hap in (1, 2):
                cache_key = (chrom, hap)
                if cache_key not in het_sites_cache:
                    try:
                        het_sites_cache[cache_key] = load_het_sites(sample, chrom, hap)
                    except FileNotFoundError:
                        # An empty dict here would NOT mean "no filter" -- het_count_in_span
                        # would return 0 for every read on every contig, silently filtering out
                        # 100% of reads at any k>=1 rather than skipping the region. Caught
                        # before this ever ran at scale: a missing-chromosome case must skip
                        # the locus entirely (matching the original unfiltered script's
                        # graceful-skip behavior for missing chain/modbed files), not proceed
                        # with a filter that looks like real data but is actually "reject
                        # everything".
                        het_sites_cache[cache_key] = None

                if het_sites_cache[cache_key] is None:
                    print(f'# SKIP {sample} hap{hap} {name} ({chrom}): no H01 het-site output '
                          f'for this (sample, chrom) -- run H01_call_hap_het_sites.py first',
                          file=sys.stderr)
                    continue

                cpg_rows, _, n_total, n_kept = parse_locus(
                    sample, hap, name, chrom, start, end,
                    het_sites_cache[cache_key], hmm_intervals, k)

                for s, h, region_name, kk, hg38_chrom, hg38_pos, n_meth, n_unmeth in cpg_rows:
                    frac = n_meth / (n_meth + n_unmeth) if (n_meth + n_unmeth) else float('nan')
                    all_rows.append((s, h, region_name, hg38_chrom, hg38_pos, n_meth, n_unmeth,
                                      frac))

                if n_total:
                    frac_kept = n_kept / n_total
                    print(f'{sample} hap{hap} {name} k={k}: {len(cpg_rows)} CpGs, '
                          f'{n_kept}/{n_total} reads kept ({frac_kept:.1%})', file=sys.stderr)

    with open(out_path, 'w') as out:
        out.write('sample\thap\tregion\thg38_chrom\thg38_pos\tn_meth\tn_unmeth\tfrac_meth\n')
        for row in all_rows:
            out.write('\t'.join(str(x) for x in row) + '\n')
    print(f'# wrote {len(all_rows)} rows to {out_path}', file=sys.stderr)


if __name__ == '__main__':
    main()
