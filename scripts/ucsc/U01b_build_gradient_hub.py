#!/usr/bin/env python3
"""U01b_build_gradient_hub.py -- second UCSC hub showing the global-methylation continuum:
hap1 per-CpG mCG for N_DONORS donors evenly spaced in hap1 global mCG (call-weighted, B01a
summaries), both chemistries (R941 + R1041; the first build used R941 only, which left a gap
between NA19338 at 0.53 and the next donor at 0.60 -- user asked to include R1041 to fill it).
R1041 runs ~4 points lower (mostly inside PMDs), so chemistry is shown in every track label.
Tracks are ordered and coloured from lowest to highest global mCG.

Hub: /u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2_gradient (same layout as U01a / the lab's hubs)
  superTrack Reference: fibroblast PMDs, LCL constitutive PMD domains (A01e), LCL PMD
                        frequency and mean mCG (10kb, copied from the asm_lr_hprc2 hub)
  superTrack Gradient:  hap1 per-CpG % mCG, one track per selected donor
Selection is written to hg38/selected_donors.tsv. Re-runnable; donors whose bigWig isn't
written yet (A02a) are reported and skipped until the next run.

Usage: U01b_build_gradient_hub.py
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
MAIN_HUB = Path('/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2/hg38')
HUB = Path('/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2_gradient')
G = HUB / 'hg38'
PER_HAP = PROJDIR / 'results' / 'meth_bins' / 'per_hap'
PMD_DIR = PROJDIR / 'data' / 'pmds'
REF_SETS = PROJDIR / 'results' / 'pmd_metagene' / 'ref_sets'
SEQ_QC = PROJDIR / 'results' / 'qc' / 'data' / 'supp_seq_qc.csv'
MANIFEST = PROJDIR / 'tsv' / 'meta' / 'hprc2_sample_manifest.tsv'
KENT = Path('/u/project/cluo/terencew/programs/kentsrc/utils')
N_DONORS = 20
N_OUTLIERS = 2
CHEMS = ('R941', 'R1041')

HUB_TXT = """hub asm_lr_hprc2_gradient
shortLabel HPRC2 LCL mCG gradient
longLabel UCLA Luo lab - HPRC2 LCL hap1 ONT CpG methylation for 20 donors (both ONT chemistries) ordered by global mCG
genomesFile genomes.txt
email cluo@mednet.ucla.edu
descriptionUrl luogenomics.github.io
"""


def viridis(x: float) -> str:
    stops = np.array([[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]], float)
    p = x * (len(stops) - 1)
    i = min(int(p), len(stops) - 2)
    c = stops[i] + (p - i) * (stops[i + 1] - stops[i])
    return ','.join(str(int(round(v))) for v in c)


def select() -> pd.DataFrame:
    chem = pd.read_csv(SEQ_QC).set_index('sample_id').sequencing_chemistry_ont
    sp = pd.read_csv(MANIFEST, sep='\t').set_index('sample_id').superpopulation
    s = pd.concat([pd.read_csv(f, sep='\t') for f in PER_HAP.glob('*_hap1.summary.tsv')])
    s = s[s['sample'].map(chem).isin(CHEMS)].sort_values('global_wmeth').reset_index(drop=True)
    # the two deepest donors (NA20762 0.415, NA19338 0.531) are far below the rest; evenly spacing
    # over the full range crowded the picks at 0.57-0.60, so keep them and space the remaining
    # picks from the 3rd-lowest donor to the top
    v = s.global_wmeth.values
    targets = np.r_[v[:N_OUTLIERS], np.linspace(v[N_OUTLIERS], v[-1], N_DONORS - N_OUTLIERS)]
    picked = []
    for t in targets:
        cand = s[~s['sample'].isin(picked)]
        picked.append(cand.loc[(cand.global_wmeth - t).abs().idxmin(), 'sample'])
    sel = s[s['sample'].isin(picked)].sort_values('global_wmeth').reset_index(drop=True)
    sel['superpopulation'] = sel['sample'].map(sp)
    sel['chemistry'] = sel['sample'].map(chem)
    sel['rank'] = np.arange(1, len(sel) + 1)
    return sel[['rank', 'sample', 'superpopulation', 'chemistry', 'global_wmeth', 'meth_in_pmd_native',
                'meth_out_pmd_native', 'mean_depth']]


def main() -> None:
    G.mkdir(parents=True, exist_ok=True)
    (HUB / 'hub.txt').write_text(HUB_TXT)
    (HUB / 'genomes.txt').write_text('genome hg38\ntrackDb hg38/trackDb.txt\n')
    shutil.copyfile(MAIN_HUB / 'hg38.chromsizes', G / 'hg38.chromsizes')
    for f in ('fibroblast_PMDs.bigBed', 'LCL_PMD_freq.bw', 'LCL_mean_mCG.bw'):
        shutil.copyfile(MAIN_HUB / f, G / f)
    tmp = G / 'lcl_constitutive.tmp.bed'
    pd.read_csv(REF_SETS / 'lcl_constitutive.bed', sep='\t', header=None).sort_values([0, 1]).to_csv(
        tmp, sep='\t', header=False, index=False)
    subprocess.run([str(KENT / 'bedToBigBed'), str(tmp), str(G / 'hg38.chromsizes'),
                    str(G / 'LCL_constitutive_PMDs.bigBed')], check=True)
    tmp.unlink()

    sel = select()
    sel.to_csv(G / 'selected_donors.tsv', sep='\t', index=False)
    print(sel.round(3).to_string(index=False))

    db = ["""track Reference
superTrack on show
shortLabel Reference
longLabel Fibroblast PMDs, LCL constitutive PMD domains, LCL PMD frequency and mean mCG (10kb)
centerLabel off
"""]
    for name, f, typ, label, color in (
            ('Fibroblast_PMDs', 'fibroblast_PMDs.bigBed', 'bigBed', 'Fibroblast PMDs (IGVF PGP start, merged, filtered)', '66,84,78'),
            ('LCL_constitutive_PMDs', 'LCL_constitutive_PMDs.bigBed', 'bigBed', 'LCL PMD domains present in >=90% of ~200 donors', '150,40,40'),
            ('LCL_PMD_freq', 'LCL_PMD_freq.bw', 'bigWig 0 1', 'Fraction of donors with this 10kb bin in their own PMD calls', '150,40,40'),
            ('LCL_mean_mCG', 'LCL_mean_mCG.bw', 'bigWig 0 1', 'Mean 10kb mCG across ~200 donors', '0,0,0')):
        extra = ('visibility dense' if typ == 'bigBed' else
                 'visibility full\n  autoScale off\n  viewLimits 0:1\n  maxHeightPixels 30:30:30\n  windowingFunction mean')
        db.append(f"""
  track {name}
  parent Reference
  type {typ}
  bigDataUrl {f}
  shortLabel {name}
  longLabel {label}
  color {color}
  {extra}
  labelFontSize 16
  longLabelFontSize 8
""")
    db.append("""
track Gradient
superTrack on show
shortLabel mCG_gradient
longLabel hap1 per-CpG mCG, 20 donors (R941+R1041) ordered from lowest (top) to highest (bottom) global mCG
centerLabel off
""")
    missing = []
    for r in sel.itertuples():
        src = PMD_DIR / r.sample / f'{r.sample}_hap1.hg38.meth.bw'
        if not src.exists():
            missing.append(r.sample)
            continue
        dst = G / f'{r.sample}_hap1_mCG.bw'
        if not dst.exists() or dst.stat().st_size != src.stat().st_size:
            shutil.copyfile(src, dst)
        color = viridis((r.rank - 1) / (N_DONORS - 1))
        db.append(f"""
  track g{r.rank:02d}_{r.sample}
  parent Gradient
  type bigWig 0 100
  bigDataUrl {dst.name}
  shortLabel {r.rank:02d}_{r.sample}
  longLabel #{r.rank} {r.sample} hap1 | {r.chemistry} | {r.superpopulation} | global mCG {r.global_wmeth:.3f} | in/out own PMDs {r.meth_in_pmd_native:.2f}/{r.meth_out_pmd_native:.2f}
  priority {r.rank}
  color {color}
  visibility full
  windowingFunction mean
  autoScale off
  yLineOnOff on
  minValue 0
  maxValue 100
  viewLimits 0:100
  maxHeightPixels 30:30:30
  labelFontSize 16
  longLabelFontSize 8
""")
    (G / 'trackDb.txt').write_text(''.join(db))
    print(f'hub written: {HUB}; missing bigWigs (rerun after A02a): {missing}')


if __name__ == '__main__':
    main()
