# asm_lr_hprc2 — population-scale allele-specific methylation from HPRC2

**Status: proposal / data-staging, 2026-09-08. Not yet an active pipeline.**

This is the successor to [`asm_lr`](https://github.com/terencewtli/asm_lr), re-platformed onto
the HPRC Release 2 (HPRC2) epigenome resource
(`epigenome.humanpangenome.org`, S3 bucket `hprc-epigenome`), which was not known to be publicly
available when `asm_lr` began. It carries forward that project's core scientific question and
its ASM-calling methodology; what changes is the input data and, as a result, what's now
statistically defensible to claim.

## Core question (unchanged from `asm_lr`)

To what extent is single-molecule methylation structure explained by genetic variation, and how
far does that genetic influence extend along the chromosome? Long reads jointly observe phased
heterozygous variants and dense CpG methylation on the same molecule, so allele-specific
methylation (ASM) can be tested at distances short-read WGBS cannot phase (~200bp limit). Each
locus is classified as **proximal genotype-driven**, **distal genotype-driven**,
**parental/imprinted**, or **non-genetic** (taxonomy from `asm_lr/md/updated_proposal.md`,
assessed against the Rosenski et al. 2025 *Nat Commun* bimodal-methylation atlas).

## What's new: why HPRC2 changes the scope

`asm_lr`'s paper scope was locked in on 2026-09-07 as biology-primary
(*"what do we gain from long-read vs. short-read ASM detection"*), ~2 months, with any
ancestry-stratified population-genetics claim **explicitly cut** — the working cohort had only
6 donors passing QC, 2 of them AFR, underpowered for any ancestry claim (the one same-ancestry
pair replicated indistinguishably from cross-ancestry pairs).

HPRC2 removes that constraint. Of 232 sample directories in `hprc-epigenome`, **229 have both
hap1 and hap2 ONT modbed files** (haplotype-resolved 5mC calls), spanning:

| Superpop | n |
|---|---|
| AFR | 63 |
| AMR | 43 |
| EAS | 49 |
| EUR | 30 |
| SAS | 36 |
| unmapped (see caveats) | 8 |

AFR — the group that was the specific statistical bottleneck in `asm_lr` (n=2) — is now the
**largest** group. A real, adequately-powered ancestry-replication analysis is possible for the
first time. This isn't scope creep: it's the population-genetics ambition the original proposal
always had, minus the data limitation that forced it out of scope one day before this resource
was found.

**All 18 of `asm_lr`'s active-cohort donors and all 30 of its original planned cohort are present
in this 229.** Nothing from the prior phase is discarded — it's a strict superset.

## What's obviated vs. what carries forward

HPRC2 hap1/hap2 modbeds are aligned directly to each donor's own diploid assembly
(contig naming `<sample>#<hap>#<contig>`, e.g. `HG00097#2#CM094075.1`) — reads are partitioned
by haplotype **because they were aligned to that haplotype's assembly**, not by downstream
statistical phasing. This means:

- **Obviated**: `asm_lr`'s W00 (align), W01/W01b (WhatsHap haplotag/phase), W02/W03 (modkit
  extract/pileup) — the entire alignment-to-pileup pipeline, including the ~38%-haplotagging-rate
  QC problem that was never resolved for the population-panel approach. HPRC2 hands us
  haplotype-partitioned calls directly.
- **Also obviated, probably**: the dipcall-vs-WhatsHap-panel phasing debate that consumed most of
  `asm_lr`'s 2026-09-07 session (`md/20260907_phasing_qc_review.md`) — assembly-based haplotype
  assignment *is* the dipcall-equivalent gold standard, provided by default for every donor.
- **Carries forward directly**: the beta-binomial ASM caller (`github/ont_asm_caller`) and its two
  confirmed calibration fixes — the region-pooling pseudoreplication correction
  (`readlevel.test_region_reads`) and the mean-dependent dispersion estimator
  (`dispersion.estimate_dispersion_trend`) — plus the validation methodology that gave the
  6-donor result its credibility: cross-donor replication as a finalization filter, checked
  against Zink et al. 2018 imprinted DMRs and Blueprint monocyte meQTLs as independent positive
  controls.

## Baseline to beat / reproduce at scale

`asm_lr`'s current 6-sane-donor result (`asm_lr/github/asm_lr/md/20260907.results.summary.md`),
to be replicated and extended at n≈229:

- 38,245 significant beta-binomial regions pooled across 6 donors → 28,255 merged loci → **3,821
  loci (13.5%) replicate in ≥2/6 donors** (finalization rule, arbitrary, not sensitivity-tested)
- **48.9%** of Zink et al. 2018 imprinted DMRs recovered by the ≥2-donor filter alone (vs. ~30%
  background replication rate)
- **~2.6×** enrichment for significant monocyte meQTLs among finalized loci with 450K probe
  coverage
- Proximal/distal split: **note the live discrepancy** — the results-summary doc above reports
  91.3%/8.2%/0.4%, but a same-day correction recorded only in project memory (not yet written
  back into that file) found this was a bug and the corrected split is closer to 50/50 (35–56%
  proximal, 42–56% distal per donor, consistent across all 6 sane donors). **Verify against
  `asm_lr` source data before quoting either number** — this is itself an instance of the
  documentation-drift problem described below.

## Scoping and data-use policy

Checked both HPRC policy pages (`humanpangenome.org/data-use/`,
`humanpangenome.org/publication-policy/`) given the concern about scooping the consortium:

- **No embargo, no pre-review, no co-authorship requirement, no "scooping" restriction.** External
  users "may freely download, analyze, and publish results... without restrictions" (AWS Open
  Data registry listing) beyond standard attribution (BioProject PRJNA730823 + NHGRI) and an
  open-access / concurrent-preprint expectation.
- **Two real constraints to respect**: no participant re-identification risk, and — directly
  relevant to a population-genetics framing — NASEM guidance to avoid unqualified continental-
  level ancestry labels. Superpopulation-level claims need careful, justified framing, not a
  bare AFR/EUR/EAS/SAS/AMR comparison presented without context.
- **The main HPRC2 preprint** (bioRxiv 2026.07.21.739710) does not perform ASM analysis — this
  project remains a genuinely orthogonal contribution, not a duplication of consortium work.

## Open items — status as of 2026-09-08, later same day

1. **modbed column format — RESOLVED.** Column layout is chrom/start/end, read ID, score,
   strand, comma-separated methylated-offset list, comma-separated unmethylated-offset list
   (Zhou et al., *Cell Genomics* 2023). The negative offsets seen in real files are **not an
   error or a quality score** — confirmed from `modbedtools` source
   (`modbed/modbed.py`, github.com/lidaof/modbedtools): the sign of each offset encodes which
   strand/stored-base the modification call came from (e.g. a call stored as `T` — the
   reverse-complement representation of a reference `C` — gets a negative offset; a call stored
   as the forward base gets a positive one). This lets one column report a CpG's calls from
   either strand at consistent genomic coordinates without a separate strand column per offset.
   Full detail in `docs/data_sources.md` §5.
2. **Assembly paths — 202/229 resolved** (was 199, now includes 3 more found under a `v1.1.0`
   version suffix instead of `v1.0.1`). **27/229 have no assembly file at all** under
   `working/HPRC/<ID>/assemblies/release2/` — not a naming mismatch, the directory is genuinely
   empty or absent. List in `docs/data_sources.md` §4.
3. **ONT basecaller/chemistry homogeneity — largely resolved, one caveat remains.** Checked raw
   ONT data folders for all 229 samples: **204/229 (89%) have a harmonized Dorado `sup5.0.0`
   re-basecall**, including **100% of the 146 samples whose original raw data was Guppy-called
   and 100% of the 56 originally Dorado-0.6-called** — i.e. the exact chemistry split that forced
   `asm_lr` to exclude 12 of its 30 donors appears to have been erased by a uniform re-basecall
   across essentially the whole cohort. **Caveat**: this confirms a shared basecaller *model* is
   available for nearly everyone, not that the epigenome bucket's `methylation.ONT.hap*.modbed.gz`
   files were themselves generated from that re-basecall rather than the older mixed one — that
   provenance link hasn't been directly confirmed. Also, a shared model can't fully erase the
   underlying R9.4.1-vs-R10.4.1 pore/signal difference, only reduce it — this de-risks pooling
   substantially but doesn't make it a non-issue to check empirically (e.g. an early QC pass
   should still look for a chemistry-batch effect in the pooled data, not assume it's gone).
   **The 25/229 samples missing the harmonized re-basecall overlap heavily with the 27 missing an
   assembly** — most likely the same underlying "not yet in this release's harmonized tier"
   stratum, and probably worth carrying as one documented exclusion set, not two.
4. **8/229 samples aren't in the 1000G 3202-sample panel** (ancestry unknown by this method) —
   cross-reference against HPRC's own sample metadata rather than assuming GIAB/non-1000G status.
5. **Sensitivity of the ≥2-donor finalization rule** was already flagged as untested in `asm_lr`
   (`md/20260907_outstanding_items.md` #11) — worth resolving before scaling, since at n≈229 the
   right replication threshold is likely not the same number tuned for n=6.

## Download infrastructure

`scripts/qsub/A01a_download_hprc2_files.sh` — SGE array job, one task per manifest row, downloads
that donor's hap1/hap2 ONT modbed (+ index) and, where resolved, hap1/hap2 assembly FASTA. See
the script header for usage; it reads `data/hprc2_sample_manifest.tsv` directly rather than
hardcoding a sample list, so re-running after the 27 unresolved/excluded samples get fixed
requires no script changes.

## Data

`data/hprc2_sample_manifest.tsv` — all 229 usable samples: population, superpopulation, ONT
modbed file sizes (hap1/hap2, as a coverage proxy), other available data types (PacBio methylC,
Fiber-seq, Hi-C, Iso-Seq — present but out of scope, see below), resolved assembly download URLs
where available, and overlap flags against `asm_lr`'s 18- and 30-donor cohorts.

`docs/data_sources.md` — how every download link and format claim above was verified (S3 listing
commands, HTTP checks run, what wasn't checked).

## Explicitly out of scope (for now)

HPRC2 also provides Hi-C, Iso-Seq expression, and Fiber-seq chromatin accessibility per donor —
exactly the "why didn't you also look at X" data that motivated the scoping concern this document
opened with. Not using them is a deliberate choice to keep this an ASM paper, not an
everything-pangenome paper; they're documented in the manifest as available in case a specific
finding later calls for one of them as a targeted follow-up, not as a standing analysis axis.
