#!/usr/bin/env python3
"""C07c_candidate_genotypes.py -- ASM replication, step 3: genotype context of every candidate in
every donor, and a per-candidate genotype-ASM association (does ASM follow heterozygosity?).

Genotypes: the cohort SNV VCF (G03, data/vcf/cohort_gw/all_donors.snps.vcf.gz), phased hap1|hap2
from the assembly-vs-hg38 calls, i.e. the same hap1/hap2 labels as the ASM delta (mu1 - mu2).
A missing genotype (./.) means the site is absent from that donor's VCF and is read as reference
(G03 merges without -0; unaligned sequence is therefore also read as reference -- a small bias
toward "homozygous" that only affects non-aligned regions).

Per variant: cohort ALT allele frequency; CpG effect vs hg38 (destroys a reference CpG: REF is
the C or G of a CG; creates one: ALT forms CG with a neighbouring base).

Per (donor, candidate):
  n_het_region      het SNVs inside the region
  n_het_5kb         het SNVs within FLANK bp of the region (incl. inside)
  nearest_het_bp    distance from the region to the nearest het SNV (0 inside; NaN if > NEAR_MAX)
  n_het_cpg_region  het SNVs inside the region that destroy or create a CpG
  lead_het, lead_delta_alt   het status at the candidate's lead variant (below) and the ASM
                    delta oriented to that variant's ALT allele (mu_ALT - mu_REF haplotype)

Per candidate (over calibrated donors -- tiers discovery + replication -- that tested it):
  cohort variant content (n SNVs / CpG-affecting SNVs in region and +-FLANK, max AF in region)
  lead variant: the SNV within +-FLANK whose heterozygosity best predicts asm_rep, one-sided
    hypergeometric p (ASM enriched in hets); needs >= 2 het and >= 2 non-het tested donors.
    Reports its AF, CpG effect, distance, the 2x2 counts, and allele-direction consistency among
    ASM hets (fraction of ASM hets whose ALT haplotype is on the majority side: ~1 for a cis
    genetic effect, ~0.5 for parent-of-origin imprinting, where hap-vs-allele is random).
  penetrance among donors with no het SNV within +-FLANK (genotype-independent ASM).

Outputs (results/asm/replication/genotype/):
  <chrom>.region.tsv.gz     one row per candidate on the chromosome
  <chrom>.donor.tsv.gz      one row per (donor, candidate) tested

Usage: C07c_candidate_genotypes.py <chrom> [outdir]   (outdir only for tests; default genotype/)
"""
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import pysam
from scipy.stats import hypergeom

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
REP = PROJDIR / 'results' / 'asm' / 'replication'
VCF = PROJDIR / 'data' / 'vcf' / 'cohort_gw' / 'all_donors.snps.vcf.gz'
HG38 = '/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosome.fa'
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'
FLANK = 5_000
NEAR_MAX = 50_000
CALIBRATED = ('discovery', 'replication')


def load_genotypes(chrom: str, windows: np.ndarray, samples: List[str]):
    # returns pos (1-based), ref, alt, hap1/hap2 ALT indicators (n_var x n_samples, int8);
    # keeps only variants inside the merged candidate windows
    cmd = [BCFTOOLS, 'query', '-r', chrom, '-f', '%POS\t%REF\t%ALT\t[%GT,]\n', str(VCF)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    ws, we = windows[:, 0], windows[:, 1]
    pos_l, ref_l, alt_l, h1_l, h2_l = [], [], [], [], []
    for chunk in pd.read_csv(proc.stdout, sep='\t', header=None, names=['pos', 'ref', 'alt', 'gt'],
                             chunksize=200_000, dtype={'pos': np.int64, 'ref': str, 'alt': str, 'gt': str}):
        p0 = chunk.pos.values - 1
        i = np.searchsorted(ws, p0, side='right') - 1
        keep = (i >= 0) & (p0 < we[np.clip(i, 0, None)])
        chunk = chunk[keep]
        if chunk.empty:
            continue
        # fixed-width 'a|b,' per sample -> chars 0 and 2 of every 4
        g = np.frombuffer(''.join(chunk['gt'].values).encode(), dtype='S1').reshape(len(chunk), -1)
        g = g[:, :4 * len(samples)].reshape(len(chunk), len(samples), 4)
        h1_l.append((g[:, :, 0] == b'1').astype(np.int8))
        h2_l.append((g[:, :, 2] == b'1').astype(np.int8))
        pos_l.append(chunk.pos.values)
        ref_l.append(chunk.ref.values)
        alt_l.append(chunk.alt.values)
    proc.wait()
    if not pos_l:
        z = np.zeros((0, len(samples)), np.int8)
        return np.zeros(0, np.int64), np.array([]), np.array([]), z, z
    return (np.concatenate(pos_l), np.concatenate(ref_l), np.concatenate(alt_l),
            np.concatenate(h1_l), np.concatenate(h2_l))


def cpg_effect(seq: str, pos1: np.ndarray, ref: np.ndarray, alt: np.ndarray) -> np.ndarray:
    # 'destroy' if REF is part of a reference CG, 'create' if ALT forms a CG, else ''
    out = np.full(len(pos1), '', dtype=object)
    for k, (p, r, a) in enumerate(zip(pos1 - 1, ref, alt)):
        prev_b = seq[p - 1] if p > 0 else 'N'
        next_b = seq[p + 1] if p + 1 < len(seq) else 'N'
        if (r == 'C' and next_b == 'G') or (r == 'G' and prev_b == 'C'):
            out[k] = 'destroy'
        elif (a == 'C' and next_b == 'G') or (a == 'G' and prev_b == 'C'):
            out[k] = 'create'
    return out


def main(chrom: str, outdir: Path = REP / 'genotype') -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    tiers = pd.read_csv(REP / 'donor_tiers.tsv', sep='\t')
    cand = pd.read_csv(REP / 'candidates.tsv.gz', sep='\t', usecols=['region_id', 'set', 'chrom', 'start', 'end'])
    cand = cand[cand.chrom == chrom].sort_values('start').reset_index(drop=True)
    print(f'{chrom}: {len(cand):,} candidates', flush=True)

    vcf_samples = subprocess.run([BCFTOOLS, 'query', '-l', str(VCF)], capture_output=True, text=True,
                                 check=True).stdout.split()
    sidx: Dict[str, int] = {s: i for i, s in enumerate(vcf_samples)}

    # merged +-NEAR_MAX windows around candidates (0-based half-open)
    w = np.stack([np.maximum(cand.start.values - NEAR_MAX, 0), cand.end.values + NEAR_MAX], axis=1)
    merged = [list(w[0])]
    for s, e in w[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    pos1, ref, alt, h1, h2 = load_genotypes(chrom, np.array(merged, dtype=np.int64), vcf_samples)
    p0 = pos1 - 1
    seq = pysam.FastaFile(HG38).fetch(chrom).upper()
    eff = cpg_effect(seq, pos1, ref, alt)
    af = (h1.sum(1) + h2.sum(1)) / (2 * h1.shape[1]) if len(pos1) else np.zeros(0)
    het = h1 != h2
    print(f'{chrom}: {len(pos1):,} SNVs in candidate windows', flush=True)

    # re-test calls for this chromosome's candidates
    cset = set(cand.region_id)
    calls = []
    for s in tiers['sample']:
        f = REP / 'calls' / f'{s}.tsv.gz'
        if f.exists():
            d = pd.read_csv(f, sep='\t', usecols=['sample', 'region_id', 'delta', 'asm_rep'])
            calls.append(d[d.region_id.isin(cset)])
    calls = pd.concat(calls, ignore_index=True).merge(tiers[['sample', 'tier']], on='sample')
    calls = calls[calls['sample'].isin(sidx)]
    by_region = {r: g for r, g in calls.groupby('region_id')}

    region_rows, donor_rows = [], []
    for r in cand.itertuples(index=False):
        lo, hi = np.searchsorted(p0, [r.start - NEAR_MAX, r.end + NEAR_MAX])
        vp = p0[lo:hi]
        dist = np.maximum(0, np.maximum(r.start - vp, vp - (r.end - 1)))
        inreg, in5 = dist == 0, dist <= FLANK
        row = {'region_id': r.region_id,
               'n_snv_region': int(inreg.sum()), 'n_snv_5kb': int(in5.sum()),
               'n_cpg_snv_region': int((eff[lo:hi][inreg] != '').sum()),
               'max_af_region': float(af[lo:hi][inreg].max()) if inreg.any() else 0.0}
        g = by_region.get(r.region_id)
        if g is None:
            region_rows.append(row)
            continue
        cols = np.array([sidx[s] for s in g['sample']])
        H = het[lo:hi][:, cols]
        n_het_region = H[inreg].sum(0)
        n_het_5kb = H[in5].sum(0)
        dmat = np.where(H, dist[:, None], np.inf)
        nearest = dmat.min(0) if len(vp) else np.full(len(cols), np.inf)
        cpgmask = inreg & (eff[lo:hi] != '')
        n_het_cpg = H[cpgmask].sum(0)

        cal = g.tier.isin(CALIBRATED).values
        asm = g.asm_rep.values.astype(bool)
        N, K = int(cal.sum()), int((asm & cal).sum())
        row.update({'n_tested_cal': N, 'n_asm_cal': K})
        nohet = (n_het_5kb == 0) & cal
        row.update({'n_nohet5kb_cal': int(nohet.sum()), 'n_nohet5kb_asm': int((nohet & asm).sum())})

        lead = -1
        if N >= 4 and K >= 1 and in5.any():
            v5 = np.flatnonzero(in5)
            Hc = H[v5][:, cal]
            n = Hc.sum(1)
            a = (Hc & asm[cal][None, :]).sum(1)
            ok = (n >= 2) & (N - n >= 2)
            if ok.any():
                p = np.where(ok, hypergeom.sf(a - 1, N, K, n), 1.0)
                j = int(np.argmin(p))
                lead = v5[j]
                row.update({'lead_pos': int(pos1[lo + lead]), 'lead_ref': ref[lo + lead], 'lead_alt': alt[lo + lead],
                            'lead_af': float(af[lo + lead]), 'lead_cpg_effect': eff[lo + lead],
                            'lead_dist': int(dist[lead]), 'lead_p': float(p[j]),
                            'lead_n_het': int(n[j]), 'lead_n_het_asm': int(a[j]),
                            'lead_n_nonhet': int(N - n[j]), 'lead_n_nonhet_asm': int(K - a[j]),
                            'n_snv_tested_5kb': int(ok.sum())})

        if lead >= 0:
            lead_het = H[lead]
            alt_on_h1 = h1[lo + lead, cols] == 1
            delta_alt = np.where(lead_het, np.where(alt_on_h1, g.delta.values, -g.delta.values), np.nan)
            m = lead_het & asm & cal
            if m.sum() >= 2:
                sgn = np.sign(delta_alt[m])
                row['lead_dir_consistency'] = float((1 + abs(sgn.mean())) / 2)
                row['lead_mean_delta_alt'] = float(np.nanmean(delta_alt[m]))
        else:
            lead_het = np.zeros(len(cols), bool)
            delta_alt = np.full(len(cols), np.nan)
        region_rows.append(row)
        donor_rows.append(pd.DataFrame({
            'sample': g['sample'].values, 'region_id': r.region_id,
            'n_het_region': n_het_region, 'n_het_5kb': n_het_5kb,
            'nearest_het_bp': np.where(np.isfinite(nearest), nearest, np.nan),
            'n_het_cpg_region': n_het_cpg, 'lead_het': lead_het, 'lead_delta_alt': delta_alt}))

    reg = cand.merge(pd.DataFrame(region_rows), on='region_id', how='left')
    for c in ('n_tested_cal', 'lead_p'):
        if c not in reg:
            reg[c] = np.nan
    reg.to_csv(outdir / f'{chrom}.region.tsv.gz', sep='\t', index=False, float_format='%.5g')
    don = pd.concat(donor_rows, ignore_index=True) if donor_rows else pd.DataFrame()
    don.to_csv(outdir / f'{chrom}.donor.tsv.gz', sep='\t', index=False, float_format='%.5g')
    tested = reg.n_tested_cal.fillna(0) > 0
    print(f'{chrom}: {tested.sum():,} candidates tested in calibrated donors; '
          f'{(reg.lead_p < 1e-3).sum():,} with lead-variant p < 1e-3; {len(don):,} donor rows', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], *(Path(a) for a in sys.argv[2:3]))
