#!/usr/bin/env python3
"""P04_read_length_coverage_qc.py -- per-sample read-length + coverage QC, for cross-donor/
cross-superpopulation/cross-PCLAI comparability checks (not ASM-specific).

Two independent sources, used for two different things -- do not conflate them:

1. Basecaller `sequencing_summary.txt.gz` (raw_data/nanopore/sup5.0.0/*_summary.txt.gz on the
   main human-pangenomics bucket) -- the STANDARD Dorado output, one row per READ with
   `sequence_length_template` and `mean_qscore_template` columns. This is the authoritative
   source for read length/quality: it's pre-alignment, pre-haplotype-split, so it isn't
   confounded by anything this project's own pipeline does. Checked first (this session) rather
   than reimplemented, matching how PCLAI/liftover chains were found instead of rebuilt.
   NOT haplotype-resolved -- a whole-donor number.

2. modbed (already downloaded, `modbed/`) -- gives aligned read SPAN per read (col3 - col2), not
   the original basecalled read length (soft-clips/indels mean these differ, usually slightly).
   Used here only for per-HAPLOTYPE coverage (sum of spans / haplotype contig length), since the
   sequencing_summary has no haplotype assignment at all.

The paper itself (2026.07.21.739710, Data section, "per-sample sequencing data and coverage are
summarized in Supplementary Table 6") implies an authoritative per-sample table already exists.
Fetching bioRxiv's supplementary-material endpoint returned HTTP 429 (rate-limited) as of
2026-09-15 -- retry that before trusting this script's numbers over the paper's own table if it
becomes reachable; this script exists as the fallback/cross-check, not a presumed replacement.

Usage: P04_read_length_coverage_qc.py <sample_id>
Prints one TSV row to stdout (header on first invocation via --header).
"""
import gzip
import sys
import urllib.request
from pathlib import Path

import numpy as np

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
MODBED_DIR = PROJDIR / 'data' / 'modbed'
CHAIN_DIR = PROJDIR / 'data' / 'chains'
BUCKET = 'https://human-pangenomics.s3.amazonaws.com'

COLUMNS = ['sample', 'n_reads_basecaller', 'median_read_length', 'read_length_n50',
           'mean_qscore', 'hap1_coverage', 'hap2_coverage']


def n50(lengths_sorted_desc, total):
    cum = 0
    for L in lengths_sorted_desc:
        cum += L
        if cum >= total / 2:
            return L
    return 0


def basecaller_read_stats(sample):
    """List raw_data/nanopore/sup5.0.0/*_summary.txt.gz for this sample via the S3 REST API
    (public, unsigned bucket -- same technique as docs/data_sources.md §1), stream each, and
    aggregate sequence_length_template + mean_qscore_template across all runs."""
    prefix = f'working/HPRC/{sample}/raw_data/nanopore/sup5.0.0/'
    list_url = f'{BUCKET}/?list-type=2&prefix={prefix}'
    with urllib.request.urlopen(list_url) as resp:
        listing = resp.read().decode()
    keys = [k for k in _extract_keys(listing) if k.endswith('_summary.txt.gz')]

    lengths = []
    qscores = []
    for key in keys:
        url = f'{BUCKET}/{key}'
        with urllib.request.urlopen(url) as resp:
            with gzip.open(resp, 'rt') as fh:
                header = fh.readline().rstrip('\n').split('\t')
                len_idx = header.index('sequence_length_template')
                q_idx = header.index('mean_qscore_template')
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    lengths.append(int(f[len_idx]))
                    qscores.append(float(f[q_idx]))

    if not lengths:
        return None
    lengths_arr = np.array(lengths)
    total = lengths_arr.sum()
    n50_val = n50(sorted(lengths, reverse=True), total)
    return dict(n_reads=len(lengths), median_length=float(np.median(lengths_arr)),
                n50=n50_val, mean_qscore=float(np.mean(qscores)))


def _extract_keys(xml_listing):
    # Minimal, dependency-free XML key extraction (same pattern used throughout this project's
    # data-discovery docs) -- avoids pulling in an XML/S3-SDK dependency for one field.
    import re
    return re.findall(r'<Key>([^<]*)</Key>', xml_listing)


def contig_sizes_from_chain(chain_path):
    """t_size per contig, read straight from chain headers (already have these files locally
    for the liftover step -- avoids a separate .fai download for a number we already have)."""
    sizes = {}
    opener = gzip.open if str(chain_path).endswith('.gz') else open
    with opener(chain_path, 'rt') as fh:
        for line in fh:
            if line.startswith('chain'):
                f = line.split()
                sizes[f[2]] = int(f[3])
    return sizes


def haplotype_coverage(sample, hap):
    modbed_path = MODBED_DIR / f'{sample}_hap{hap}.modbed.gz'
    chain_path = CHAIN_DIR / f'{sample}_hap{hap}_vs_GRCh38.chain.gz'
    if not modbed_path.exists() or not chain_path.exists():
        return None
    sizes = contig_sizes_from_chain(chain_path)
    genome_size = sum(sizes.values())
    total_aligned_bases = 0
    with gzip.open(modbed_path, 'rt') as fh:
        for line in fh:
            f = line.split('\t', 3)
            total_aligned_bases += int(f[2]) - int(f[1])
    return total_aligned_bases / genome_size if genome_size else None


def main():
    if sys.argv[1:2] == ['--header']:
        print('\t'.join(COLUMNS))
        return
    sample = sys.argv[1]
    stats = basecaller_read_stats(sample)
    hap1_cov = haplotype_coverage(sample, 1)
    hap2_cov = haplotype_coverage(sample, 2)
    row = [sample,
           stats['n_reads'] if stats else '',
           f"{stats['median_length']:.1f}" if stats else '',
           stats['n50'] if stats else '',
           f"{stats['mean_qscore']:.2f}" if stats else '',
           f'{hap1_cov:.2f}' if hap1_cov is not None else '',
           f'{hap2_cov:.2f}' if hap2_cov is not None else '']
    print('\t'.join(str(x) for x in row))


if __name__ == '__main__':
    main()
