#!/usr/bin/env python3
"""C08d_imprinted_domains.py -- reclassify candidates against imprinted DOMAINS, not a 10 kb
window around individual Zink DMRs.

Why. C07d flags imprinting as `dist_zink <= 10 kb`. Imprinted loci are domains, not points:
SNRPN/PWS spans ~2 Mb, KCNQ1OT1 ~1 Mb, ZDBF2/GPR1-AS ~300 kb. The 10 kb rule therefore splits a
single imprinted domain into an `imprinting` core plus a ring of regions that fall through to
`genotype_indep`. Measured: of the 155 candidates at penetrance >= 0.5 outside the 10 kb window,
127 are classed `genotype_indep`, their nearest genes are ZDBF2, PWAR1, SNORD116-30, SNHG14 and
H19, and their lead-variant direction consistency is 0.579 -- i.e. ~0.5, parent-of-origin, not
cis. They are imprinting.

Domains are built by merging Zink DMRs within MERGE_GAP of each other and padding by PAD, then
(diagnostically) intersected with a curated imprinted-gene name list so the domain calls can be
sanity-checked by eye rather than trusted blind.

Direction consistency is also promoted to a first-class discriminator, because it is a MECHANISM
test rather than an overlap test and, unlike `genotype_linked`'s 1.00, it is not circular for
these classes: parent-of-origin gives ~0.5, cis gives ~1.

Outputs (results/asm/model/):
  imprinted_domains.bed            the merged domains
  candidates_reclassified.tsv.gz   region_id + rep_class (C07d) + rep_class_v2 + domain id
  reclass_summary.tsv              old x new contingency and per-class features

Usage: C08d_imprinted_domains.py [--merge-gap 1000000] [--pad 100000]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
REP = PROJDIR / 'results' / 'asm' / 'replication'
REGIONS = PROJDIR / 'results' / 'asm' / 'regions' / 'regions_hg38.tsv.gz'
OUT = PROJDIR / 'results' / 'asm' / 'model'

# canonical imprinted loci, for labelling the merged domains (not for defining them)
IMPRINTED = ['SNRPN', 'SNURF', 'SNHG14', 'PWAR1', 'PWRN1', 'PWRN2', 'PWRN4', 'UBE3A', 'MAGEL2',
             'NDN', 'MKRN3', 'IPW', 'H19', 'IGF2', 'KCNQ1', 'KCNQ1OT1', 'CDKN1C', 'PHLDA2',
             'GNAS', 'GNASAS1', 'NESP', 'MEST', 'MESTIT1', 'GRB10', 'PEG3', 'PEG10', 'ZIM2',
             'ZDBF2', 'GPR1', 'GPR1-AS', 'PLAGL1', 'HYMAI', 'DLK1', 'MEG3', 'MEG8', 'RTL1',
             'DIO3', 'NNAT', 'BLCAP', 'L3MBTL1', 'NAP1L5', 'INPP5F', 'ZNF331', 'ZNF597',
             'NAA60', 'ERLIN2', 'MCTS2', 'HM13', 'FAM50B', 'TFPI2', 'SGCE', 'PPP1R9A',
             'CALCR', 'DDC', 'HERC3', 'NAP1L4', 'OSBPL5', 'SLC22A18', 'TSSC4', 'TRPM5',
             'CPA4', 'THSD7A', 'GLIS3', 'WT1', 'ANO1', 'RB1', 'AIRN', 'IGF2R', 'SLC22A2',
             'SLC22A3', 'PLAG1', 'NTM', 'DIRAS3', 'TP73', 'LRRTM1', 'COPG2', 'KLF14',
             'SNORD108', 'SNORD109A', 'SNORD115', 'SNORD116', 'MIR298']


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--merge-gap', type=int, default=1_000_000)
    ap.add_argument('--pad', type=int, default=100_000)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    reg = pd.read_csv(REGIONS, sep='\t')
    z = reg[reg['set'] == 'zink'][['chrom', 'start', 'end']].sort_values(['chrom', 'start'])
    print('%d Zink DMRs -> merging within %d bp, padding %d bp' % (len(z), a.merge_gap, a.pad), flush=True)

    doms = []
    for chrom, g in z.groupby('chrom'):
        s = g.start.values - a.pad
        e = g.end.values + a.pad
        cur_s, cur_e, k = s[0], e[0], 1
        for i in range(1, len(s)):
            if s[i] <= cur_e + a.merge_gap:
                cur_e = max(cur_e, e[i])
                k += 1
            else:
                doms.append((chrom, max(cur_s, 0), cur_e, k))
                cur_s, cur_e, k = s[i], e[i], 1
        doms.append((chrom, max(cur_s, 0), cur_e, k))
    dom = pd.DataFrame(doms, columns=['chrom', 'start', 'end', 'n_dmr'])
    dom['domain_id'] = ['imprdom_%s_%d' % (c, s) for c, s in zip(dom.chrom, dom.start)]
    dom['span_kb'] = (dom.end - dom.start) / 1000.0
    dom.to_csv(OUT / 'imprinted_domains.bed', sep='\t', index=False, header=False,
               columns=['chrom', 'start', 'end', 'domain_id', 'n_dmr'])
    print('%d imprinted domains; span median %.0f kb, max %.0f kb; %d multi-DMR' % (
        len(dom), dom.span_kb.median(), dom.span_kb.max(), (dom.n_dmr > 1).sum()), flush=True)

    c = pd.read_csv(REP / 'candidates_replication.tsv.gz', sep='\t', low_memory=False)
    c['in_impr_domain'] = False
    c['domain_id'] = ''
    for chrom, g in dom.groupby('chrom'):
        m = c.chrom == chrom
        if not m.any():
            continue
        mid = ((c.loc[m, 'start'] + c.loc[m, 'end']) // 2).values
        idx = np.searchsorted(g.end.values, mid, side='left')
        idx = np.clip(idx, 0, len(g) - 1)
        inside = (mid >= g.start.values[idx]) & (mid < g.end.values[idx])
        c.loc[m, 'in_impr_domain'] = inside
        c.loc[m, 'domain_id'] = np.where(inside, g.domain_id.values[idx], '')

    dirc = c.lead_dir_consistency
    pofo = dirc <= 0.65                      # ~0.5: parent-of-origin
    cis = dirc >= 0.8                        # ~1.0: cis genetic
    replicating = c.rep_class.isin(['imprinting', 'genotype_linked', 'genotype_indep', 'other'])

    v2 = c.rep_class.astype(object).copy()
    # promote: replicating + inside an imprinted domain + parent-of-origin direction -> imprinting
    promote = replicating & c.in_impr_domain & (pofo | c.zink_10kb)
    v2[promote] = 'imprinting'
    # replicating, in a domain, but clearly cis -> keep the genetic label
    v2[replicating & c.in_impr_domain & cis & ~c.zink_10kb] = 'genotype_linked'
    c['rep_class_v2'] = v2

    c[['region_id', 'chrom', 'start', 'end', 'rep_class', 'rep_class_v2', 'in_impr_domain',
       'domain_id', 'penetrance_rep', 'lead_dir_consistency', 'zink_10kb', 'nearest_gene']] \
        .to_csv(OUT / 'candidates_reclassified.tsv.gz', sep='\t', index=False, float_format='%.5g')

    ct = pd.crosstab(c.rep_class, c.rep_class_v2)
    ct.to_csv(OUT / 'reclass_summary.tsv', sep='\t')
    print('\n=== old (rows) x new (cols) ===')
    print(ct.to_string())
    hi = c[c.penetrance_rep >= 0.5]
    print('\n=== regions at penetrance >= 0.5 (n=%d) ===' % len(hi))
    print('  old class:', hi.rep_class.value_counts().to_dict())
    print('  new class:', hi.rep_class_v2.value_counts().to_dict())
    print('  frac imprinting: %.3f -> %.3f' % (
        (hi.rep_class == 'imprinting').mean(), (hi.rep_class_v2 == 'imprinting').mean()))
    moved = c[(c.rep_class != c.rep_class_v2)]
    print('\n%d regions reclassified; their median penetrance %.3f, median dir_consistency %.3f' % (
        len(moved), moved.penetrance_rep.median(), moved.lead_dir_consistency.median()))
    print('  top nearest genes among reclassified:',
          moved.nearest_gene.value_counts().head(10).to_dict())
    known = moved.nearest_gene.isin(IMPRINTED)
    print('  %.1f%% of reclassified regions sit next to a canonically imprinted gene' % (100 * known.mean()))
    print('\nwrote', OUT / 'candidates_reclassified.tsv.gz')


if __name__ == '__main__':
    main()
