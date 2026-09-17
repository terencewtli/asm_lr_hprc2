#!/usr/bin/env python3
"""U01a_build_track_hub.py -- UCSC track hub for asm_lr_hprc2 at
/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2, same layout as the lab's other hubs there
(hub.txt, genomes.txt, hg38/trackDb.txt with superTracks, real files, no symlinks).

Tracks (hg38):
  LCL_population  -- fibroblast PMDs (igvf_pgp, copied from ../fibro_d0), and 10kb maps across
                     all donors: mean mCG, SD of mCG, PMD frequency, caller-free low-methylation
                     frequency, PMD-boundary frequency (B01b / B02a)
  <donor> x2      -- NA19338 (lowest global mCG) and HG04187 (highest): hg38 PMDs per haplotype
                     and hap-consensus (B01a), per-CpG mCG per haplotype (A02a, 0-100)
  LCL_donors_10kb -- every donor haplotype's 10kb mean mCG (B01a, 0-1), hidden by default;
                     blue = R941, orange = R1041 (HPRC2 Supp S6); labels carry superpopulation

Per-CpG bigWigs are ~720 MB each, so only the two highlighted donors get them (all 404 would be
~290 GB). Re-runnable: rebuilds everything from current results; per-CpG tracks are added once
A02a has written them.

Usage: U01a_build_track_hub.py
"""
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2')
HUB = Path('/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2')
G = HUB / 'hg38'
BINS = PROJDIR / 'results' / 'meth_bins'
PER_HAP = BINS / 'per_hap'
PMD_DIR = PROJDIR / 'data' / 'pmds'
FIB_BIGBED = Path('/u/project/cluo/PUBLIC_SHARED/ucsc/fibro_d0/hg38/start_merged.filt_mcg.sorted.bigBed')
FULL_SIZES = Path('/u/project/cluo/PUBLIC_SHARED/ucsc/fibro_d0/hg38/hg38.chromsizes')
AUTO_SIZES = Path('/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.autosomal.chromsizes')
SEQ_QC = PROJDIR / 'results' / 'qc' / 'data' / 'supp_seq_qc.csv'
MANIFEST = PROJDIR / 'tsv' / 'meta' / 'hprc2_sample_manifest.tsv'
KENT = Path('/u/project/cluo/terencew/programs/kentsrc/utils')
BEDTOOLS = '/u/home/t/terencew/bin/bedtools'
HIGHLIGHT = {'NA19338': 'lowest global mCG in cohort; flagged by HPRC2 methylation QC',
             'HG04187': 'highest global mCG in cohort'}
CHEM_COLOR = {'R941': '40,81,176', 'R1041': '230,120,20'}
BIN = 10_000

HUB_TXT = """hub asm_lr_hprc2
shortLabel HPRC2 LCL ONT mCG
longLabel UCLA Luo lab - HPRC2 lymphoblastoid haplotype-resolved ONT CpG methylation, PMDs and population domain maps
genomesFile genomes.txt
email cluo@mednet.ucla.edu
descriptionUrl luogenomics.github.io
"""


def sizes() -> pd.DataFrame:
    s = pd.read_csv(AUTO_SIZES, sep='\t', header=None, names=['chrom', 'len'])
    return s[s.chrom.isin([f'chr{i}' for i in range(1, 23)])]


def bins_to_bw(df: pd.DataFrame, col: str, out: Path) -> None:
    sz = sizes()
    ln = dict(zip(sz.chrom, sz.len))
    bw = pyBigWig.open(str(out), 'w')
    bw.addHeader(list(zip(sz.chrom, sz.len.astype(int))))
    df = df.dropna(subset=[col])
    for c in sz.chrom:
        g = df[df.chrom == c].sort_values('bin_start')
        if not len(g):
            continue
        s = g.bin_start.values.astype(np.int64)
        e = np.minimum(s + BIN, ln[c])
        bw.addEntries([c] * len(g), s.tolist(), ends=e.tolist(), values=g[col].astype(float).tolist())
    bw.close()


def bed_to_bigbed(bed: pd.DataFrame, out: Path) -> bool:
    if not len(bed):
        return False
    tmp = out.with_suffix('.tmp.bed')
    bed.sort_values(['chrom', 'start']).to_csv(tmp, sep='\t', header=False, index=False)
    subprocess.run([str(KENT / 'bedToBigBed'), str(tmp), str(FULL_SIZES), str(out)], check=True)
    tmp.unlink()
    return True


def read_pmd(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=['chrom', 'start', 'end'])
    return pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])


def bw_stanza(name: str, parent: str, url: str, label: str, color: str, vis: str, lo: float,
              hi: float, height: int = 40) -> str:
    return f"""
  track {name}
  parent {parent}
  type bigWig
  bigDataUrl {url}
  shortLabel {name}
  longLabel {label}
  color {color}
  visibility {vis}
  windowingFunction mean
  autoScale off
  yLineOnOff on
  minValue {lo}
  maxValue {hi}
  viewLimits {lo}:{hi}
  maxHeightPixels {height}:{height}:{height}
  labelFontSize 16
  longLabelFontSize 8
"""


def bb_stanza(name: str, parent: str, url: str, label: str, color: str) -> str:
    return f"""
  track {name}
  parent {parent}
  type bigBed
  bigDataUrl {url}
  shortLabel {name}
  longLabel {label}
  color {color}
  visibility dense
  labelFontSize 16
  longLabelFontSize 8
"""


def super_stanza(name: str, label: str) -> str:
    return f"""
track {name}
superTrack on show
shortLabel {name}
longLabel {label}
autoScale on
aggregate transparentOverlay
centerLabel off
"""


def main() -> None:
    G.mkdir(parents=True, exist_ok=True)
    (HUB / 'hub.txt').write_text(HUB_TXT)
    (HUB / 'genomes.txt').write_text('genome hg38\ntrackDb hg38/trackDb.txt\n')
    shutil.copyfile(FULL_SIZES, G / 'hg38.chromsizes')
    shutil.copyfile(FIB_BIGBED, G / 'fibroblast_PMDs.bigBed')
    db = []

    # population maps
    db.append(super_stanza('LCL_population', 'HPRC2 LCL population maps (10kb bins, ~200 donors) and fibroblast PMDs'))
    db.append(bb_stanza('Fibroblast_PMDs', 'LCL_population', 'fibroblast_PMDs.bigBed',
                        'Fibroblast PMDs (IGVF PGP start, merged across four donors, filtered)', '66,84,78'))
    freq = pd.read_csv(BINS / 'domain_frequency_10kb.tsv.gz', sep='\t')
    for col, name, label, color, hi in (
            ('mean_meth', 'LCL_mean_mCG', 'Mean 10kb mCG across donors (donor = mean of haps)', '0,0,0', 1),
            ('sd_meth', 'LCL_sd_mCG', 'SD of 10kb mCG across donors', '120,120,120', 0.2),
            ('freq_pmd', 'LCL_PMD_freq', 'Fraction of donors with this 10kb bin in their own genome-wide PMD calls', '150,40,40', 1),
            ('freq_rel', 'LCL_lowmCG_freq', 'Fraction of donors with bin mCG < donor median - 0.10 (caller-free)', '200,90,40', 1)):
        bins_to_bw(freq, col, G / f'{name}.bw')
        db.append(bw_stanza(name, 'LCL_population', f'{name}.bw', label, color, 'full', 0, hi))
    bf = pd.read_csv(BINS / 'boundaries' / 'boundary_freq_10kb.tsv.gz', sep='\t')
    bins_to_bw(bf, 'freq', G / 'LCL_PMD_boundary_freq.bw')
    db.append(bw_stanza('LCL_PMD_boundary_freq', 'LCL_population', 'LCL_PMD_boundary_freq.bw',
                        'Fraction of donors with a hap-consensus PMD boundary within +/-10kb', '90,40,150', 'full', 0, 1))

    # highlighted donors
    for d, why in HIGHLIGHT.items():
        db.append(super_stanza(d, f'{d} ({why}): PMDs and per-CpG mCG per haplotype'))
        beds = {h: read_pmd(PER_HAP / f'{d}_hap{h}.hg38.pmd.bed') for h in (1, 2)}
        for h, bed in beds.items():
            if bed_to_bigbed(bed, G / f'{d}_hap{h}_PMDs.bigBed'):
                db.append(bb_stanza(f'{d}_hap{h}_PMDs', d, f'{d}_hap{h}_PMDs.bigBed',
                                    f'{d} hap{h} PMDs (genome-wide dnmtools pmd, native coords, CpG-mapped to hg38)', '150,40,40'))
        if len(beds[1]) and len(beds[2]):
            a, b = (G / f'{d}_tmp{h}.bed' for h in (1, 2))
            for h, p in ((1, a), (2, b)):
                beds[h].sort_values(['chrom', 'start']).to_csv(p, sep='\t', header=False, index=False)
            inter = subprocess.run([BEDTOOLS, 'intersect', '-a', str(a), '-b', str(b)], check=True,
                                   capture_output=True, text=True).stdout
            a.unlink()
            b.unlink()
            cons = pd.DataFrame([l.split('\t')[:3] for l in inter.strip().split('\n') if l],
                                columns=['chrom', 'start', 'end']).astype({'start': int, 'end': int})
            if bed_to_bigbed(cons, G / f'{d}_consensus_PMDs.bigBed'):
                db.append(bb_stanza(f'{d}_consensus_PMDs', d, f'{d}_consensus_PMDs.bigBed',
                                    f'{d} PMDs present on both haplotypes', '100,0,0'))
        for h in (1, 2):
            src = PMD_DIR / d / f'{d}_hap{h}.hg38.meth.bw'
            if src.exists():
                shutil.copyfile(src, G / f'{d}_hap{h}_mCG.bw')
                db.append(bw_stanza(f'{d}_hap{h}_mCG', d, f'{d}_hap{h}_mCG.bw',
                                    f'{d} hap{h} per-CpG mCG (%, CpG-filtered ONT calls, all reads)',
                                    '40,81,176' if h == 1 else '176,40,81', 'full', 0, 100))
            else:
                print(f'# {src.name} not written yet; rerun after A02a finishes')

    # all donors, 10kb
    chem = pd.read_csv(SEQ_QC).set_index('sample_id').sequencing_chemistry_ont
    sp = pd.read_csv(MANIFEST, sep='\t').set_index('sample_id').superpopulation
    summ = pd.concat([pd.read_csv(f, sep='\t') for f in PER_HAP.glob('*.summary.tsv')]).set_index(['sample', 'hap'])
    db.append(super_stanza('LCL_donors_10kb', 'Per-haplotype 10kb mean mCG, every donor (blue R941, orange R1041)'))
    (G / 'donors_10kb').mkdir(exist_ok=True)
    for f in sorted(PER_HAP.glob('*.bins10kb.tsv.gz')):
        tag = f.name[:-len('.bins10kb.tsv.gz')]
        s, h = tag.rsplit('_hap', 1)
        out = G / 'donors_10kb' / f'{tag}.bw'
        if not out.exists() or out.stat().st_mtime < f.stat().st_mtime:
            bins_to_bw(pd.read_csv(f, sep='\t', usecols=['chrom', 'bin_start', 'mean_meth']), 'mean_meth', out)
        c = chem.get(s, 'NA')
        gm = summ.loc[(s, int(h)), 'global_wmeth'] if (s, int(h)) in summ.index else np.nan
        vis = 'full' if s in HIGHLIGHT else 'hide'
        db.append(bw_stanza(tag, 'LCL_donors_10kb', f'donors_10kb/{tag}.bw',
                            f'{s} hap{h} 10kb mCG | {sp.get(s, "NA")} | {c} | global mCG {gm:.3f}',
                            CHEM_COLOR.get(c, '100,100,100'), vis, 0, 1, 25))
    (G / 'trackDb.txt').write_text(''.join(db))
    n = sum(1 for _ in G.glob('donors_10kb/*.bw'))
    print(f'hub written: {HUB} ({n} donor-hap 10kb tracks)')


if __name__ == '__main__':
    main()
