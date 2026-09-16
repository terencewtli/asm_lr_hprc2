#!/usr/bin/env python3
"""H04_epiallele_classification.py -- disambiguates what H03's flat purity-vs-k result
actually means, per (region, hap): is a "low purity" locus genuinely bimodal (two real,
separated read populations -- misassignment or true epiallele mixture is a coherent
explanation there) or is it uniform/intermediate (every read individually partially
methylated, no real two-state structure at all -- purity's binary >=0.5 threshold is just
the wrong tool there, not evidence of a phasing problem)?

Same GMM(k=1) vs GMM(k=2) BIC-comparison idea as asm_lr's existing
notebooks/old/chr19/A07_epiallele_complexity.ipynb, applied here to this project's own
per-read frac_meth values (from H02's *.reads.tsv) instead of re-deriving read-level fractions
from scratch -- reuses the concept, not the chr19-specific code (that notebook is tied to a
different data layout: pickled epiallele_store per chr19 locus, not this project's modbed-
derived tables).

Usage: H04_epiallele_classification.py <reads_tsv> <out_tsv>
  reads_tsv: the *.reads.tsv file H02_filtered_locus_matrix.py writes (needs >= MIN_READS
    reads per region x hap to classify at all)
"""
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

MIN_READS = 8  # GMM(k=2) needs more data to be a meaningful fit than the simple majority-vote
               # purity metric did (MIN_READS=3 there) -- too few points and a 2-component fit
               # is just overfitting noise, not detecting real structure.


def classify(frac_meth_values):
    x = np.array(frac_meth_values).reshape(-1, 1)
    gmm1 = GaussianMixture(n_components=1, random_state=0).fit(x)
    gmm2 = GaussianMixture(n_components=2, random_state=0).fit(x)
    bic1, bic2 = gmm1.bic(x), gmm2.bic(x)
    is_bimodal = bic2 < bic1
    # for a bimodal fit, report how separated the two components are (in frac_meth units) --
    # a "bimodal" fit with near-identical component means isn't really two states
    if is_bimodal:
        means = sorted(gmm2.means_.flatten())
        separation = means[1] - means[0]
    else:
        separation = 0.0
    return dict(bic1=bic1, bic2=bic2, is_bimodal=is_bimodal, separation=separation,
                mean=float(np.mean(x)), std=float(np.std(x)))


def main():
    reads_tsv, out_tsv = sys.argv[1], sys.argv[2]
    reads = pd.read_csv(reads_tsv, sep='\t')

    rows = []
    for (region, hap), g in reads.groupby(['region', 'hap']):
        vals = g['frac_meth'].values
        if len(vals) < MIN_READS:
            continue
        result = classify(vals)
        n_high = (vals >= 0.5).sum()
        purity = max(n_high, len(vals) - n_high) / len(vals)
        rows.append(dict(region=region, hap=hap, n_reads=len(vals), purity=purity, **result))

    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep='\t', index=False)

    n_bimodal = df['is_bimodal'].sum()
    print(f'{len(df)} (region, hap) pairs classified (>= {MIN_READS} reads)', flush=True)
    print(f'{n_bimodal} ({n_bimodal/len(df):.1%}) genuinely bimodal (GMM k=2 preferred by BIC)',
          flush=True)
    print(f'{len(df)-n_bimodal} ({1-n_bimodal/len(df):.1%}) unimodal/intermediate -- purity\'s '
          f'"disagreement" there is not real two-state structure, the metric is measuring the '
          f'wrong thing', flush=True)
    print(flush=True)
    print('mean purity, bimodal vs unimodal loci:', flush=True)
    print(df.groupby('is_bimodal')['purity'].agg(['mean', 'count']).to_string(), flush=True)
    print(flush=True)
    print('among bimodal loci, component separation (0=identical, 1=fully separated):',
          flush=True)
    print(df[df['is_bimodal']]['separation'].describe().to_string(), flush=True)


if __name__ == '__main__':
    main()
