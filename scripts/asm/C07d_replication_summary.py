#!/usr/bin/env python3
"""C07d_replication_summary.py -- ASM replication, step 4: penetrance, replication classes, and
the feature comparison of replicating vs non-replicating candidates.

Penetrance is computed on the REPLICATION tier only (1.2 <= lambda_gc <= 3, genomic-controlled,
candidate-set BH; C07b) -- those donors had no part in choosing the candidates. Also reported:
the same within R9.4.1 donors only (chemistry check, RESULTS §1 rule), and discovery-tier counts
(circular, for reference only).

Classes (cpg candidates):
  untested       n_tested_rep < MIN_TESTED
  private        no replication donor calls it
  sporadic       called in some replication donors, but not more than a background rate P0 would
                 give (binomial p >= P_REP)
  replicating    binomial p < P_REP against P0, then split by mechanism:
    imprinting        within 10 kb of a Zink imprinted DMR
    genotype_linked   lead-variant het/ASM association p < P_LEAD and allele-direction
                      consistency >= DIR_MIN (a cis genetic effect)
    genotype_indep    penetrance among donors with no het SNV within 5 kb >= 0.5 of overall
                      penetrance (ASM without nearby sequence difference: imprinting-like or
                      epigenetic, not yet annotated)
    other             replicating, none of the above
P0 is the median per-donor asm_rep rate over all candidates in replication donors, i.e. what a
donor calls at a random candidate -- deliberately conservative (it includes true signal).

Outputs (results/asm/replication/):
  candidates_replication.tsv.gz   candidates + penetrance + genotype summary (C07c) + class
  class_summary.tsv               per class: n and feature summaries (imprinting, mQTL, CGI-like,
                                  CpG o/e, TSS distance, PMD class, RT, lead AF, CpG-SNV lead)
  chromhmm_by_class.tsv           ChromHMM group composition per class
  nearest_het_by_call.tsv         nearest-het distance distribution for ASM vs non-ASM donor
                                  calls, per class (replication donors)

Usage: C07d_replication_summary.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
REP = PROJDIR / 'results' / 'asm' / 'replication'
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
MIN_TESTED = 10
P_REP = 1e-3
P_LEAD = 1e-4
DIR_MIN = 0.8


def main() -> None:
    tiers = pd.read_csv(REP / 'donor_tiers.tsv', sep='\t')
    cand = pd.read_csv(REP / 'candidates.tsv.gz', sep='\t')
    rep_donors = tiers[tiers.tier == 'replication']
    r941 = set(rep_donors[rep_donors.chemistry == 'R941']['sample'])

    calls = pd.concat([pd.read_csv(REP / 'calls' / f'{s}.tsv.gz', sep='\t',
                                   usecols=['sample', 'region_id', 'asm_rep'])
                       for s in tiers['sample'] if (REP / 'calls' / f'{s}.tsv.gz').exists()])
    calls = calls.merge(tiers[['sample', 'tier']], on='sample')
    rc = calls[calls.tier == 'replication']
    p0 = rc.groupby('sample').asm_rep.mean().median()
    print(f'replication donors with calls: {rc["sample"].nunique()}; background P0 = {p0:.4f}', flush=True)

    pen = rc.groupby('region_id').asm_rep.agg(n_tested_rep='size', n_asm_rep='sum')
    pen941 = rc[rc['sample'].isin(r941)].groupby('region_id').asm_rep.agg(n_tested_rep941='size', n_asm_rep941='sum')
    pdisc = calls[calls.tier == 'discovery'].groupby('region_id').asm_rep.agg(n_tested_disc='size', n_asm_disc='sum')
    d = cand.merge(pen, left_on='region_id', right_index=True, how='left') \
            .merge(pen941, left_on='region_id', right_index=True, how='left') \
            .merge(pdisc, left_on='region_id', right_index=True, how='left')
    for c in ['n_tested_rep', 'n_asm_rep', 'n_tested_rep941', 'n_asm_rep941', 'n_tested_disc', 'n_asm_disc']:
        d[c] = d[c].fillna(0).astype(int)
    with np.errstate(all='ignore'):
        d['penetrance_rep'] = d.n_asm_rep / d.n_tested_rep
        d['penetrance_rep941'] = d.n_asm_rep941 / d.n_tested_rep941
    d['p_replicate'] = binom.sf(d.n_asm_rep - 1, d.n_tested_rep, p0)

    geno = pd.concat([pd.read_csv(REP / 'genotype' / f'{c}.region.tsv.gz', sep='\t')
                      for c in AUTOSOMES if (REP / 'genotype' / f'{c}.region.tsv.gz').exists()])
    geno = geno.drop(columns=['set', 'chrom', 'start', 'end'])
    d = d.merge(geno, on='region_id', how='left')
    with np.errstate(all='ignore'):
        pen_nohet = d.n_nohet5kb_asm / d.n_nohet5kb_cal
        pen_cal = d.n_asm_cal / d.n_tested_cal
    d['penetrance_nohet5kb'] = pen_nohet

    rep = d.p_replicate < P_REP
    cls = np.select(
        [d.n_tested_rep < MIN_TESTED, d.n_asm_rep == 0, ~rep,
         d.zink_10kb,
         (d.lead_p < P_LEAD) & (d.lead_dir_consistency >= DIR_MIN),
         (d.n_nohet5kb_cal >= 3) & (pen_nohet >= 0.5 * pen_cal)],
        ['untested', 'private', 'sporadic', 'imprinting', 'genotype_linked', 'genotype_indep'],
        default='other')
    d['rep_class'] = np.where(d.set == 'zink', 'zink_control', cls)
    d.to_csv(REP / 'candidates_replication.tsv.gz', sep='\t', index=False, float_format='%.5g')

    order = ['zink_control', 'imprinting', 'genotype_linked', 'genotype_indep', 'other',
             'sporadic', 'private', 'untested']
    g = d.groupby('rep_class')
    summ = pd.DataFrame({
        'n': g.size(),
        'median_n_tested_rep': g.n_tested_rep.median(),
        'median_penetrance_rep': g.penetrance_rep.median(),
        'median_penetrance_rep941': g.penetrance_rep941.median(),
        'frac_zink_10kb': g.zink_10kb.mean(),
        'frac_mqtl_overlap': g.mqtl_overlap.mean(),
        'frac_cgi_like': g.cgi_like.mean(),
        'median_cpg_obs_exp': g.cpg_obs_exp.median(),
        'median_cpg_per_100bp': g.cpg_per_100bp.median(),
        'median_dist_tss': g.dist_tss.median(),
        'frac_tss_2kb': g.dist_tss.apply(lambda x: (x <= 2000).mean()),
        'frac_constitutive_pmd': g.class_pmd.apply(lambda x: (x == 'constitutive').mean()),
        'median_rt': g.rt.median(),
        'median_snv_5kb': g.n_snv_5kb.median(),
        'frac_cpg_snv_region': g.n_cpg_snv_region.apply(lambda x: (x > 0).mean()),
        'median_lead_af': g.lead_af.median(),
        'frac_lead_cpg_effect': g.lead_cpg_effect.apply(lambda x: x.isin(['destroy', 'create']).mean()),
        'median_lead_dist': g.lead_dist.median(),
        'median_lead_dir_consistency': g.lead_dir_consistency.median(),
    }).reindex([o for o in order if o in g.groups])
    summ.to_csv(REP / 'class_summary.tsv', sep='\t', float_format='%.4g')
    ch = pd.crosstab(d.rep_class, d.chromhmm_group, normalize='index').reindex(summ.index)
    ch.to_csv(REP / 'chromhmm_by_class.tsv', sep='\t', float_format='%.4g')

    # nearest het SNV for ASM vs non-ASM calls, replication donors, per class
    dons = []
    for c in AUTOSOMES:
        f = REP / 'genotype' / f'{c}.donor.tsv.gz'
        if f.exists():
            dons.append(pd.read_csv(f, sep='\t', usecols=['sample', 'region_id', 'nearest_het_bp', 'n_het_5kb']))
    don = pd.concat(dons).merge(rc, on=['sample', 'region_id']).merge(d[['region_id', 'rep_class']], on='region_id')
    bins = [-1, 0, 500, 1000, 5000, 20000, 50000, np.inf]
    labels = ['inside', '1-500', '500-1k', '1-5k', '5-20k', '20-50k', '>50k/none']
    don['nearest_het_bin'] = pd.cut(don.nearest_het_bp.fillna(np.inf), bins, labels=labels)
    nh = pd.crosstab([don.rep_class, don.asm_rep], don.nearest_het_bin, normalize='index')
    nh.to_csv(REP / 'nearest_het_by_call.tsv', sep='\t', float_format='%.4g')

    pd.set_option('display.width', 250)
    print('\nclass summary:\n' + summ.T.round(3).to_string())
    print('\nnearest het SNV (rows: class, asm call):\n' + nh.round(3).to_string())


if __name__ == '__main__':
    main()
