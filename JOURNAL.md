# Journal

Distilled project state. **Read this first.** The full session-by-session record (including every
bug, wrong turn and retraction) is in `JOURNAL.archive.md`, date-ordered, verbatim as of
2026-09-16. This file keeps only:

1. **Status board** — what's been checked and can be trusted, what's retracted, what's open.
2. **Planned analyses.**
3. **Log** — short dated entries (newest first). Put long narrative in the archive, not here.

Labels: **[verified]** = rechecked directly against data, with the check described;
**[reported]** = produced by an earlier session and not independently rechecked since;
**[retracted]** = shown wrong, kept so nobody reuses it.

---

## 1. Status board (as of 2026-09-17)

### Data provenance and quality

- **[verified] Modbeds come from the WashU HPRC Epigenome Browser bucket**
  (`s3://hprc-epigenome/samples/<ID>/methylation.ONT.hap{1,2}.modbed.gz`, via
  `scripts/download/A01a_download_hprc2_files.sh`). They are read-level: one row per read, CpG
  offsets in cols 7/8, already split into hap1/hap2 by aligning reads to each haplotype assembly.
  The mapping of offsets to positions is correct: ~90% of calls land on a CG in the assembly, for
  both strands.
- **[verified] ONT data is NOT uniformly basecalled, and the modbeds use the ORIGINAL calls.**
  This corrects the 2026-09-08 "chemistry largely resolved" conclusion.
  - Every donor's `raw_data/nanopore/` folder falls into one of these groups:

    | original data | donors | extra `sup5.0.0*` folder |
    |---|---|---|
    | R9.4.1 + Guppy 6.x, 5mC only | 146 | new R10.4.1 runs (sequenced 2024) |
    | R10.4.1 + Dorado 0.6 / sup4.3, 5mC + 5hmC | 56 | a re-basecall of the same runs |
    | no ONT BAMs in the bucket (HPP partner sites?) | 25 | none |
    | only `sup5.0.0` R10 | 2 | — |

    So for R9 donors, "has sup5" means newer R10 data exists, not that the R9 data was harmonized.
  - NA19338 and HG00099: 6,000/6,000 sampled modbed read IDs each match the donor's 2022
    **R9.4.1 Guppy** runs; none match its 2024 R10 sup5 runs. Provenance for R10 donors can't be
    settled by read ID, because the dorado0.6 and sup5 folders contain the same reads.
- **[verified] HPRC2 documents chemistry per sample:** Supp Table S6 `sequencing_chemistry_ont`
  lists 156 R941 and 73 R1041. It agrees with the folder-based labels for all 202 classifiable
  donors, and each donor's released coverage comes from one chemistry. **Use S6 as the canonical
  chemistry variable** (copied to `results/qc/data/supp_seq_qc.csv`).
  - The preprint's promoter-mQTL model includes "ONT chemistry" as a covariate, alongside 30 PEER
    factors, 10 genetic PCs, coverage and N50. So the authors know about it.
  - The main text reports no chemistry effect size and no global-methylation analysis.
  - Their analysis repo (github.com/twlab/HPRC2_DNA_Methylation) returns 404 as of 2026-09-17.
  - The user plans to contact the WashU authors.
- **[verified] Chemistry is a major global-methylation covariate, bigger than ancestry**
  (219 donors, QC04 global methylation):
  - Variance explained: R² = 0.219 for chemistry alone (R9 vs R10 vs other), 0.116 for
    superpopulation alone, 0.265 for both.
  - R10 donors are about 3.2 points lower (0.623 vs 0.659). The gap appears within every
    superpopulation (AFR 0.613 vs 0.655, SAS 0.631 vs 0.668, …).
  - After adjusting for chemistry, superpopulation effects are ≤1.5 points (AMR, SAS, p ≈ 0.03–0.05).
  - Chemistry is confounded with superpopulation: EUR is 19/27 R10, AMR 36/41 R9.
  - 13 of the 20 lowest-methylation donors are R10.
  - One candidate mechanism, not tested: R10/Dorado calls 5hmC separately, while R9 Guppy
    5mC-only calls may absorb 5hmC.
- **[verified] Low-methylation donors appear in every superpopulation.** Donor counts with global
  methylation < 0.60: AFR 10, EUR 4, EAS 4, SAS 2. Minimum per superpopulation: AFR 0.523,
  EUR 0.558, EAS 0.580, SAS 0.587, AMR 0.621. Any per-superpopulation summary must be adjusted
  for chemistry first.
- **[verified] HPRC2's own QC flags NA19338, NA20762 and HG02583** as having consistently abnormal
  ONT methylation across most or all of their BAMs (preprint Methods, "Nanopore Methylation QC").
  NA19338 is our reference PMD donor, and NA20762 is another PMD-positive donor.
- **[verified] Per-CpG depth is ~30x per haplotype** (NA19338 hap1: mean/median 30.1/30 over
  32.2M assembly CpGs, 99.93% of all assembly CpGs; ~35x physical read depth, mean read span
  57kb). Coverage is the same in R9 and R10 donors (median 28 vs 29 on chr20).
- **[verified] LCL provenance** (preprint Methods + Supp Table S15):
  - All lines are from Coriell: 210 established there, 23 established elsewhere.
  - Expanded to 4×10⁸ cells, i.e. passage 5 for 170 lines; the rest are p3–p11, 28 have no record.
  - Three ONT extraction/library protocols were used; 14 samples came from HPP partner sites.
  - No culture timeframe is given.
- **[reported]** 202/229 donors have a release-2 assembly; HG00272 has no chain file;
  8 samples are unmapped to a superpopulation.

### Methcounts / PMD pipeline

- **[retracted] "Median per-CpG coverage = 1; ONT per-CpG call depth is much sparser than read depth"**
  (archive, 2026-09-16 "still later").
  - Cause: `A01a_modbed_to_methcounts.py` wrote every called position, including non-CpG junk
    positions at depth ~1, which outnumber real CpGs 2:1 (67.8M junk vs 32.2M real rows).
  - Fixed by writing only positions where the assembly reads `CG`; output is now
    `*.cpg.methcounts.tsv.gz`.
  - The diploid-pooling plan built on that claim is dropped.
- **[retracted] "All 202 donors processed."** The arrays only ran manifest rows 1–202, which
  covered 177 of the donors with assemblies. They now run all 458 tasks and skip donors with
  missing inputs.
- **[verified] The junk-row fix barely changes PMD calls.** NA19338 chr20 hap1: 72 PMDs covering
  63.6% of chr20 (vs 64.2% before), old-vs-new Jaccard 0.946.
- **[verified] Genome-wide `dnmtools pmd` (A01c) calls PMDs in EVERY haplotype**:
  1.24–1.72 Gb per haplotype, median 1.45 Gb (≈40–56% of the assembly), domains up to ~15 Mb
  (first 81 haplotypes; HG00097 hap1: 40% of CpGs in PMDs, mean methylation 0.56 inside vs 0.73
  outside).
  - Fit on a whole genome, the HMM always finds a lower and a higher compartment: a relative call.
  - The chr20-only fits below behave differently (most donors get 0), so **chr20-only PMD burden
    numbers are not comparable to genome-wide ones.**
  - Consequence: use continuous per-bin methylation (B01) as the primary readout, and treat PMD
    calls as a relative compartment label.
- **[verified] chr20 PMD survey, 201 donors × 2 haps** (`results/qc/data/chr20_pmd_survey/`;
  PMDs called on the chr20 contig alone):
  - Each haplotype uses its single best chr20 contig; a few donors have chr20 split across
    contigs (~310K instead of ~730K CpGs).
  - **Only 8 donors get any chr20-only `dnmtools pmd` calls, and calls are all-or-nothing**:
    51–65% of chr20 or 0%.
    - The 8: NA19338 (R9), HG03017 (R9), HG00128, HG03784, NA19682, NA20282, NA20762, NA21093
      (all R10). That's 6/51 R10 donors vs 2/142 R9 donors.
    - HG03017 and NA21093 get calls on hap2 only.
  - **Every donor (197/197) is less methylated inside NA19338's PMDs than outside** (difference
    0.06–0.36), including the 10 most-methylated donors (0.06–0.10), all of which get 0 calls.
  - **Methylation inside the PMD regions is continuous and tracks global methylation**:
    r = 0.995 with global methylation (188 donors); the inside/outside gap widens as global
    methylation falls (r = −0.95).
    - The calls don't follow this: HG00290, HG03521 and HG02583 are more hypomethylated in
      these regions (0.43–0.45) than called donor HG00128, yet get 0 calls.
    - **Conclusion: the domains are present cohort-wide; `dnmtools pmd` thresholds a continuum.**
      PMD call burden is not a usable per-donor phenotype.
    - After adjusting for global methylation, chemistry has no further effect on in-PMD
      methylation (p = 0.28).
  - **Haplotypes:** mean |hap1 − hap2| in 10kb bins is 0.026 inside the reference PMDs vs 0.019
    outside. The haplotypes are nearly equal, consistent with the user's short-read WGBS result.
  - NA19338 chr20 PMDs: mean methylation 0.403 inside vs 0.653 outside; median size 349kb,
    mean 578kb, largest 2.70Mb.
- **[verified] LCL PMD locations are highly reproducible across donors; LCL vs fibroblast is not**
  (chr20, hg38, consensus-of-haps; `boundaries/`):
  - LCL vs LCL Jaccard is 0.78–0.92 for the 5 donors with full calls (both R9 NA19338 and the R10
    donors). HG03784, with fewer calls, is 0.40–0.44.
  - LCL vs fibroblast (`reference/igvf_pgp`) Jaccard is 0.41–0.43 (0.25 for HG03784). The
    fibroblast set is largely nested inside the LCL set (79–91% of fibroblast PMD bp are inside
    LCL PMDs), and LCL PMDs cover about twice as much sequence.
  - Boundaries: the median distance from an LCL boundary to the nearest fibroblast boundary is
    119–233kb. 55/126 of NA19338's boundaries are within 10kb of a fibroblast boundary.
  - Gene density at boundary windows (±5kb) is about 1.0 genes vs 0.82 for random 10kb windows.
    That's weak and untested on a single chromosome.
  - Genes at boundaries shared with fibroblasts: CD40, SIRPB2, NKX2-4, XRN2, SRSF6, …
  - Genes at LCL-specific boundaries seen in ≥4 donors: JAG1, EYA2, SULF2, KCNB1, SPO11, …
  - Genes at fibroblast-specific boundaries: BMP2, KCNQ2, SNAP25, PLCB4, MACROD2, SLC24A3, …
  - chr20 only, no statistics; wait for the genome-wide run before interpreting.
- **[retracted/superseded]** QC06's threshold-based PMD percentages. The QC10 fibroblast-overlap
  numbers came from pilot calls on junk-inflated input; their chr20 values reproduce (Jaccard
  0.41), but they should be recomputed genome-wide.
- **[reported]** `dnmtools pmd -S/-r/-p` error on this build; only `-o` is used.
- **[verified] Vectorized chain mapping replaces liftOver for per-CpG lifting** (B01a).
  - Speed: liftOver needed ~3h per haplotype for ~32M single-base records; the vectorized
    version takes 3 min and 8.8 GB (HG00097 hap1).
  - Agreement on NA19338 chr20: reproduces 99.3% of liftOver positions with identical values,
    and 98.4% land on an hg38 CG (liftOver: 97.6%). The difference is the strand-flipped-block
    correction (the native C maps to the hg38 G, so the CpG is at mapped − 1).
  - A02a bigWigs and the chr20 survey used liftOver, which misses that correction, so CpGs on
    flipped blocks are 1 bp off there. Flipped blocks carry <1% of chain score; this is minor.
  - Genome-wide mapping rate is 81.8% (chr20: 92%); centromere/satellite sequence doesn't map.

### HPRC2 preprint (bioRxiv 2026.07.21.739710) — what it does and doesn't cover [verified, main text + Methods + supp xlsx]
- **Methylation content:** the "A Panepigenome" section covers:
  - 17.6M non-reference CpGs;
  - PacBio vs ONT agreement (1.55% mean deviation over 1kb bins; details in Supp Note 2, not
    local);
  - graph-based methylation calling (Panmethyl);
  - 80,854 promoter mQTLs (Supp S11), with local-ancestry allele-frequency differences of lead
    variants (PM20D1 example);
  - a Nanopore methylation QC that flags NA20762, NA19338 and HG02583.
- **Not in the main text:**
  - PMDs or large hypomethylated domains (no "partially methylated" hits);
  - global-methylation variance decomposition by ancestry, sex, age or passage;
  - a chemistry effect size (chemistry appears only as an mQTL covariate);
  - an LCL clonality discussion.
  - Supplementary Notes/Figs weren't available locally, so this is not an exhaustive check.
- **Supp tables:**
  - S10: CpG counts and methylation by genomic context and variant source.
  - S11: *significant* promoter mQTLs only. 200bp T2T-CHM13 bins with lead variant and q; no
    effect sizes, no tested-but-null bins, no genome-wide var-CpG table.
  - S15: passage and karyotype.
  - S6: chemistry and coverage.
- **S11 lifted to hg38:** 80,016/80,854 bins, using UCSC `hs1ToHg38.over.chain.gz`, stored in
  `/u/project/cluo_scratch/terencew/claude/asm_lr_hprc2/ref/`.
- **chr20 pilot:** 52.7% of mQTL bins fall inside NA19338's consensus PMDs, vs 44% of gene-end
  ±1kb sequence (61% of chr20 overall). This hints that genetic regulation isn't depleted in
  PMDs, but the tested-bin background isn't published.

### ASM / phasing pipeline

- **[reported]** H01–H04 (2026-09-15):
  - Het sites come from assembly-vs-assembly alignment.
  - The `k≥1` het-site filter is the default; raising k gives no purity gain (≤0.08).
  - Only 32.7% of chr15 imprinted-DMR loci are bimodal.
- **[reported]** Het-site density is strongly ancestry-structured (R² = 0.917), so raw
  cross-ancestry ASM detection rates are confounded. Documented in `ont_asm_caller` PRIORITIES
  item 8.
- **[reported]** Post-reorg path bugs were fixed (2026-09-16). P03 production matrix job 14764741
  is still running.
- **[reported]** Dipcall-style VCFs (G01→G03, jobs 14766207/10/13): validated on 5 donors
  (ts/tv 1.94); full-run completion not rechecked.
- **[reported]** QC09 chain-gap asymmetry: mean 0.6% of the genome (max 1.43%).
- **[reported]** Ancestry PCA: PC1 51%, PC2 18.8%. Output is in `reference/1000G/pca/asm_lr_hprc2/`.

### Running jobs (submitted 2026-09-16)

| job | what | depends on |
|---|---|---|
| 14772522 | A01a CpG-only methcounts, all contigs, 458 tasks | — |
| 14772523 | A01c genome-wide `dnmtools pmd` (native haplotype coordinates) | 14772522 |
| 14772524 | A02a per-haplotype hg38 bigWigs (per-CpG % methylation, autosomes) | 14772522 |
| 14764741 | P03 production donor×chrom ASM matrix | — |
| 14773897 | B01a genome-wide hg38 10kb bins + hg38 PMD intervals per haplotype (`scripts/meth_bins/`) | A01c |
| 14774195 | B01b_summary: B01b matrix/PCA/covariates/continuum/domain-frequency + B02a genome-wide boundary genes | B01a |

The chr20 survey (14772521) and A01a (14772522) are done: 404/404 CpG methcounts. The user can now
delete the old unfiltered methcounts (154G):
`rm /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*_hap?.methcounts.tsv.gz`

### Open decisions

1. **Chemistry handling.** Per S6: 156 R941 and 73 R1041.
   - Re-basecalling and re-mapping is ruled out: no bandwidth in a ~2-month project (user,
     2026-09-17). The user will ask the WashU authors how they handle it.
   - Plan: S6 chemistry is a covariate in every cross-donor model; every headline result is
     replicated within R941 alone (n≈150) and shown per chemistry.
   - For ASM specifically: hap1 and hap2 of a donor share chemistry, so the within-donor ASM
     contrast is largely self-controlled. Chemistry still matters through:
     - per-donor noise and dispersion (estimate the caller's nuisance parameters per donor,
       which already happens, and report them by chemistry);
     - 5mC-vs-5hmC separation (R10 Dorado splits 5hmC out, R9 Guppy 5mC may absorb it), which
       can shift absolute haplotype methylation at 5hmC-rich loci;
     - any cross-donor "penetrance" or frequency statistic, which needs chemistry adjustment and
       an R941-only replication.
   - Positive controls to check: imprinted-DMR ASM detection rate and effect size by chemistry.
2. **When to lift to hg38.** Currently PMDs are called in native coordinates and only per-CpG
   values are lifted. For cross-donor work, lifting CpGs first gives common coordinates and a
   common CpG set, at the cost of ~8% of CpGs that don't lift. Suggest doing both.
3. **Framing (recommendation, 2026-09-17; user to decide).**
   - Don't make this Figure 1 of the ASM paper as a PMD story. Make it a Figure 1 panel set
     titled "what the ONT haplotype methylomes look like and what varies across donors",
     because every downstream ASM/penetrance claim depends on it:
     - chemistry (S6) as the dominant technical axis;
     - a continuum of global methylation / domain hypomethylation that is shared across donors
       at the same locations and is not ancestry-structured after chemistry adjustment;
     - near-equal hap1/hap2 at domains (so domain state is not a source of false ASM);
     - HPRC2-flagged donors.
   - The genome-wide PCA, domain-frequency map and LCL-vs-fibroblast comparison go to
     supplementary figures, plus a Methods/QC note.
   - If domain frequency shows constitutive LCL domains plus a donor-state continuum, and
     genetics-vs-state variance partitions cleanly by domain class (planned analysis D), that is
     an LCL-methylome result HPRC2 didn't report. It could become a short separate note or
     resource paper rather than stretching the ASM paper's theme.
   - Decide after B01b/B02a land.
4. **Superpopulation summaries.** Present them only after adjusting for chemistry; lead with
   donor-level distributions.
5. The 25 donors with no ONT BAMs in the bucket have unknown chemistry; the 27 donors without an
   assembly are excluded.

---

## 2. Planned analyses

**A. Sample × genomic-bin methylation matrix.**
- Windows of 10kb (and 50kb) × ~200 donors, per haplotype and pooled.
- Per bin: mean CpG methylation, number of CpGs, coverage, and variance / fraction of reads
  methylated.
- Inputs: the A02a hg38 bigWigs or lifted CpG methcounts.
- Analysis: PCA and hierarchical clustering.
- Question: what are the main axes of variation (global methylation / PMD state, ancestry,
  chemistry batch, something else)? Regress PCs on chemistry, superpopulation, global
  methylation, passage and HPRC2 QC flags.

**B. Is there a global methylation continuum?**
- Per donor: global CpG methylation, fraction of the genome with low methylation, PMD burden,
  median methylation in predefined PMD-like regions, and the number and total length of
  low-methylation domains. Plot them against each other.
- The chr20 survey already suggests a continuum: PMD-region methylation tracks global
  methylation at r = 0.995, while the HMM calls are discrete.
- Next: repeat genome-wide with a threshold-free domain measure, e.g. PMD-region vs
  flanking-region methylation, or a per-bin low-methylation fraction.

**C. Are domain locations stable across donors?** (the most interesting one)
- For each 10kb hg38 window, compute the fraction of donors in which it is hypomethylated
  (relative to that donor's own distribution, so it isn't just a global shift). This gives a
  population frequency map.
- Classify windows as:
  - **constitutive**: 90–100% of LCLs;
  - **variable**: 20–80%;
  - **individual-specific**: a few donors.
- This tells whether NA19338 is (A) an extreme version of a common LCL methylome or (B) a
  qualitatively unusual state. chr20 so far points to (A): 0.78–0.92 Jaccard among PMD-called
  donors, and every donor is lower in the same regions.
- Then compare the constitutive/variable sets and their boundaries to the fibroblast
  reprogramming PMDs (igvf_pgp) and to gene annotations. Boundary genes are expected to be
  enriched for interesting genes; test genome-wide against matched random windows.

**A–C are implemented genome-wide** in `scripts/meth_bins/`:
- B01a: per haplotype.
- B01b: matrix, PCA (raw and per-donor-centered, 10kb and 50kb), PC~covariate R², continuum
  metrics, domain frequency (own-PMD and caller-free relative definitions, per chemistry),
  donor×donor and donor×fibroblast PMD-bin Jaccard.
- B02a: genome-wide LCL vs fibroblast boundary loci, genes, and permutation enrichment.

**D. Genetically regulated methylation vs PMDs** (worth doing; design):
- For HPRC2 S11 promoter mQTL bins (lifted to hg38): compare the rate inside constitutive,
  variable and never-PMD 10kb bins against a matched background of all promoter 200bp bins
  (built from gencode TSSs; the tested-bin universe isn't published, so match on CpG density
  and gene-TSS proximity).
- With our own data: per-bin between-donor variance decomposed into genetic (het-site/ASM-based,
  via P03) vs global-state (regression on donor global methylation) components, by domain
  class. Predicted: PMD variance is dominated by donor state, CGI/promoter variance by genetics.

**Also pending:**
- Genome-wide version of the chr20 boundary/gene analysis (A03d) once 14772523 lands.
- Recompute QC10.
- Chemistry-stratified replication of every per-donor result above.

---

## 3. Log (newest first)

### 2026-09-17
- Genome-wide A01c PMD calls cover ~40–56% of every haplotype, so the HMM is relative.
  chr20-only burden numbers are not comparable to genome-wide ones.
- Found S6 `sequencing_chemistry_ont` (156 R941 / 73 R1041) and confirmed the preprint uses ONT
  chemistry as an mQTL covariate. Reviewed the preprint for PMDs, variance decomposition and
  var-CpGs (none reported; S11 has significant promoter mQTLs only, in CHM13). Lifted S11 to
  hg38 and ran the chr20 mQTL-in-PMD pilot.
- Built `scripts/meth_bins/`:
  - B01a: per haplotype; vectorized chain mapping validated against liftOver, ~3 min/hap.
  - B01b: matrix, PCA, covariates, continuum, domain frequency, bin-level Jaccards.
  - B02a: genome-wide boundary genes.
  - All tested on 4 donors (outputs will be overwritten by the full run).
  - Jobs 14773897 (B01a array, held on A01c) and 14774195 (summary, held on B01a).
- chr20 survey caveats noted: single best contig; liftOver misses the flipped-block offset.
- Added ASM chemistry-handling plan, planned analysis D (mQTL/genetic vs state variance by
  domain class), and the framing recommendation.

### 2026-09-16 (late)
- Retracted the median-1 coverage claim (the junk-row bug). Added the CpG filter; fixed the
  array range (458 tasks) and three A02a bigWig crashes. Resubmitted A01a, A01c and A02a.
- chr20 cross-donor survey (A03a/b) → PMD regions form a cohort-wide continuum, the caller is
  all-or-nothing, and haplotypes are nearly equal.
- hg38 PMD intervals from CpG runs (A03c) and the LCL-vs-fibroblast boundary/gene comparison
  (A03d).
- Found that the modbeds use original R9/R10 calls (not the sup5 harmonization) and that
  chemistry explains R² = 0.22 of global methylation. HPRC2's QC flags NA19338.
- Journal split into this distilled file plus `JOURNAL.archive.md`.
- Scripts: `scripts/call_pmds/{all_donors,chr20_survey}/`, `scripts/bigwig/`. Mirrored to GitHub
  (`sync_to_github.sh` now covers both).

### ≤ 2026-09-16 (earlier sessions)
See `JOURNAL.archive.md`:
- 2026-09-08: project founded, downloads.
- 2026-09-15: haplotype-assignment / het-site filter.
- 2026-09-16: path-reorg fixes, popgen QC sweep, PMD pilot, HMM PMD port, and the now-retracted
  coverage finding.

**Entry template (keep entries short; long narrative goes to the archive):**
```
### YYYY-MM-DD — summary
- what changed / what was verified (with the check)
- jobs, files
- status-board lines added/changed above
```
