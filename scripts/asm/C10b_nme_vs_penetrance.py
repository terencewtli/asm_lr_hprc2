#!/usr/bin/env python3
"""C10b_nme_vs_penetrance.py -- the one pre-registered test C10a exists to run.

PRE-REGISTERED (written before the C10a array finished, JOURNAL 2026-09-19):

  H1  Epiallele entropy predicts penetrance. Among `genotype_linked` loci, per-locus mean NME
      correlates with `pen_het` (penetrance among lead-variant heterozygotes). Direction not
      fixed a priori -- a SWITCH-like element (low entropy, clean two-state) could go either way
      against a GRADED one (high entropy, reads individually intermediate) -- so this is a
      two-sided test and the direction is part of the result, not part of the hypothesis.

  H2  Replication of `asm_lr`'s W12 pilot (n=1 donor): ASM near a het SNV is MORE deterministic,
      i.e. lower within-haplotype entropy, than ASM far from one. Directional, one-sided.

  H3  Within a locus, donors who CALL ASM differ in NME from donors who do not. This is the
      donor-resolved version of H1 and is the one that could show entropy is a state rather than
      a locus property.

Calibration gate, checked FIRST and reported whatever it says: `t_nme` must exceed `nme_null`
(hap1's reads split by whole molecule and compared to itself). `pattern.py` has a simulation
benchmark but no real-data validation, so if the observed statistic does not separate from its
own null, none of H1-H3 are interpretable and the section is not written.

Outputs (results/asm/model/):
  nme_per_locus.tsv.gz        per-locus NME summary joined to penetrance and class
  nme_vs_penetrance.tsv       the three tests, effect sizes, and the calibration gate

Usage: C10b_nme_vs_penetrance.py
"""
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
REP = PROJDIR / 'results' / 'asm' / 'replication'
MODEL = PROJDIR / 'results' / 'asm' / 'model'


def main() -> None:
    fs = sorted(glob.glob(str(MODEL / 'epiallele' / '*.nme.tsv.gz')))
    print('donors with NME: %d' % len(fs), flush=True)
    d = pd.concat([pd.read_csv(f, sep='\t') for f in fs], ignore_index=True)
    print('donor x locus rows: %d ; loci: %d' % (len(d), d.region_id.nunique()), flush=True)

    res = []

    # ---- calibration gate ----
    obs, nul = d.t_nme.dropna(), d.nme_null.dropna()
    u = stats.mannwhitneyu(obs, nul, alternative='greater')
    ratio = obs.median() / max(nul.median(), 1e-12)
    print('\n=== CALIBRATION GATE ===')
    print('  t_nme median %.4f vs null median %.4f  -> %.2fx  (MWU p=%.3g)' % (
        obs.median(), nul.median(), ratio, u.pvalue))
    gate = (ratio > 1.2) and (u.pvalue < 1e-6)
    print('  gate %s' % ('PASSED' if gate else 'FAILED -- H1-H3 are not interpretable'))
    res.append({'test': 'calibration_gate', 'stat': ratio, 'p': u.pvalue, 'pass': gate,
                'note': 't_nme median / nme_null median'})

    # ---- per-locus summary ----
    per = d.groupby('region_id').agg(
        n_donors=('sample', 'size'), nme_mean=('nme1', 'mean'), nme2_mean=('nme2', 'mean'),
        t_nme=('t_nme', 'median'), nme_null=('nme_null', 'median'),
        n_cpg=('n_cpg', 'median'), rep_class=('rep_class', 'first')).reset_index()
    per['nme_locus'] = (per.nme_mean + per.nme2_mean) / 2

    c7 = pd.read_csv(REP / 'candidates_replication.tsv.gz', sep='\t',
                     usecols=['region_id', 'lead_n_het', 'lead_n_het_asm', 'penetrance_rep',
                              'lead_dist', 'cpg_obs_exp', 'rep_class'], low_memory=False)
    c7['pen_het'] = c7.lead_n_het_asm / c7.lead_n_het
    per = per.merge(c7.drop(columns=['rep_class']), on='region_id', how='left')
    per.to_csv(MODEL / 'nme_per_locus.tsv.gz', sep='\t', index=False, float_format='%.5g')

    # ---- H1 ----
    gl = per[(per.rep_class == 'genotype_linked') & per.pen_het.notna() & per.nme_locus.notna()]
    r1 = stats.spearmanr(gl.nme_locus, gl.pen_het)
    print('\n=== H1: NME vs penetrance among genotype_linked (n=%d) ===' % len(gl))
    print('  spearman = %+.3f  p = %.3g' % (r1.correlation, r1.pvalue))
    print('  (for scale, §13 reports the best structural predictor at rho = +0.235)')
    gl2 = gl.copy()
    gl2['q'] = pd.qcut(gl2.nme_locus, 4, labels=['Q1 low entropy', 'Q2', 'Q3', 'Q4 high entropy'])
    print(gl2.groupby('q').agg(n=('region_id', 'size'), med_pen_het=('pen_het', 'median')).round(3).to_string())
    res.append({'test': 'H1_nme_vs_pen_het', 'stat': r1.correlation, 'p': r1.pvalue,
                'pass': r1.pvalue < 0.05, 'note': 'genotype_linked, n=%d' % len(gl)})

    # ---- H2 ----
    near = per[(per.lead_dist <= 500) & per.nme_locus.notna()].nme_locus
    far = per[(per.lead_dist > 5000) & per.nme_locus.notna()].nme_locus
    if len(near) > 20 and len(far) > 20:
        u2 = stats.mannwhitneyu(near, far, alternative='less')
        print('\n=== H2: asm_lr W12 replication -- is proximal ASM more deterministic? ===')
        print('  NME near lead (<=500bp, n=%d) %.4f  vs  far (>5kb, n=%d) %.4f ; one-sided p=%.3g'
              % (len(near), near.median(), len(far), far.median(), u2.pvalue))
        print('  %s' % ('REPLICATES (proximal lower entropy)' if u2.pvalue < 0.05
                        else 'does NOT replicate at n=%d donors' % len(fs)))
        res.append({'test': 'H2_proximal_determinism', 'stat': float(near.median() - far.median()),
                    'p': u2.pvalue, 'pass': u2.pvalue < 0.05, 'note': 'near<=500bp vs far>5kb'})

    # ---- H3 ----
    calls = pd.concat([pd.read_csv(REP / 'calls' / ('%s.tsv.gz' % Path(f).name.split('.')[0]),
                                   sep='\t', usecols=['sample', 'region_id', 'asm_rep'])
                       for f in fs if (REP / 'calls' / ('%s.tsv.gz' % Path(f).name.split('.')[0])).exists()],
                      ignore_index=True)
    dd = d.merge(calls, on=['sample', 'region_id'], how='inner')
    dd['nme_mean'] = (dd.nme1 + dd.nme2) / 2
    if dd.asm_rep.nunique() == 2:
        # within-locus contrast: centre NME on the locus mean so this is not the H1 effect again
        dd['nme_c'] = dd.nme_mean - dd.groupby('region_id').nme_mean.transform('mean')
        a = dd.loc[dd.asm_rep == 1, 'nme_c'].dropna()
        b = dd.loc[dd.asm_rep == 0, 'nme_c'].dropna()
        u3 = stats.mannwhitneyu(a, b)
        print('\n=== H3: within-locus, ASM-calling vs non-calling donors ===')
        print('  locus-centred NME: ASM+ %+.4f (n=%d) vs ASM- %+.4f (n=%d) ; p=%.3g' % (
            a.median(), len(a), b.median(), len(b), u3.pvalue))
        res.append({'test': 'H3_within_locus_asm_vs_not', 'stat': float(a.median() - b.median()),
                    'p': u3.pvalue, 'pass': u3.pvalue < 0.05, 'note': 'locus-centred NME'})

    pd.DataFrame(res).to_csv(MODEL / 'nme_vs_penetrance.tsv', sep='\t', index=False, float_format='%.5g')
    print('\nwrote', MODEL / 'nme_vs_penetrance.tsv')


if __name__ == '__main__':
    main()
