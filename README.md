# asm_lr_hprc2 — population-scale allele-specific methylation from HPRC2

**Status: proposal / data-staging, 2026-09-08. Not yet an active pipeline.**

Successor to [`asm_lr`](https://github.com/terencewtli/asm_lr), re-platformed onto the HPRC
Release 2 (HPRC2) epigenome resource (`epigenome.humanpangenome.org`, S3 bucket
`hprc-epigenome`). Same scientific question and ASM-calling methodology as `asm_lr`; what
changes is the input data and what's now statistically defensible to claim. Full history of why
this project exists and the scope-lock it followed → `docs/20260908_project_transition.md`.

## Core question (unchanged from `asm_lr`)

To what extent is single-molecule methylation structure explained by genetic variation, and how
far does that genetic influence extend along the chromosome? Long reads jointly observe phased
heterozygous variants and dense CpG methylation on the same molecule, enabling ASM tests at
distances short-read WGBS can't phase (~200bp limit). Each locus is classified as **proximal
genotype-driven**, **distal genotype-driven**, **parental/imprinted**, or **non-genetic**
(taxonomy from `asm_lr/md/updated_proposal.md`, assessed against Rosenski et al. 2025 *Nat
Commun* bimodal-methylation atlas).

## Why HPRC2

`asm_lr` locked its scope to biology-primary (long-read vs. short-read ASM detection) on
2026-09-07, explicitly cutting any ancestry-stratified claim — the working cohort had only 6 QC-
passing donors, 2 of them AFR, underpowered. HPRC2 removes that constraint: 229/232 sample
directories have both hap1 and hap2 ONT modbed files, spanning:

| Superpop | n |
|---|---|
| AFR | 63 |
| AMR | 43 |
| EAS | 49 |
| EUR | 30 |
| SAS | 36 |
| unmapped | 8 |

AFR (the old bottleneck, n=2) is now the **largest** group — a real, adequately-powered ancestry-
replication analysis is possible for the first time. All 18 of `asm_lr`'s active-cohort donors
and all 30 of its original planned cohort are present in this 229 (strict superset — nothing
prior is discarded).

## Pipeline: what's obviated, what carries forward

HPRC2 hap1/hap2 modbeds are aligned directly to each donor's own diploid assembly — haplotype
partitioning comes from assembly alignment, not downstream statistical phasing.

- **Obviated**: `asm_lr`'s align → WhatsHap haplotag/phase → modkit extract/pileup pipeline
  (W00–W03), including its ~38%-haplotagging-rate QC problem.
- **Carries forward**: the beta-binomial ASM caller (`github/ont_asm_caller`) with its two
  confirmed calibration fixes (region-pooling pseudoreplication correction, mean-dependent
  dispersion estimator), plus the cross-donor-replication validation methodology (checked
  against Zink et al. 2018 imprinted DMRs and Blueprint monocyte meQTLs as independent positive
  controls).

`asm_lr`'s 6-donor baseline to reproduce/extend at n≈229, and the current open QC/format items
for HPRC2, are tracked in `docs/20260908_project_transition.md`.

## Scoping and data-use

No embargo, no pre-review, no co-authorship requirement per HPRC's data-use/publication policy
— standard attribution only. One real constraint for a population-genetics framing: NASEM
guidance against unqualified continental-level ancestry labels, so superpopulation-level claims
need careful, justified framing. Hi-C, Iso-Seq, Fiber-seq, and PacBio methylC are also available
per donor but out of scope — this stays an ASM paper, not an everything-pangenome paper.

## Data & infra

- `data/hprc2_sample_manifest.tsv` — all 229 usable samples: population/superpopulation, ONT
  modbed file sizes (coverage proxy), other available data types, resolved assembly download
  URLs where available, overlap flags vs. `asm_lr`'s 18- and 30-donor cohorts.
- `docs/data_sources.md` — how every download link and format claim was verified.
- `scripts/qsub/A01a_download_hprc2_files.sh` — SGE array download job; reads the manifest
  directly rather than hardcoding a sample list.
- `docs/20260908_project_transition.md` — full detail on the scope-lock history, the baseline
  numbers to beat, open QC/format items, and known documentation-drift caveats as of 2026-09-08.
