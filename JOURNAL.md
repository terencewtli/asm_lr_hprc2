# Journal

Distilled project state. **Read this first.** The full session-by-session record (including every
bug, wrong turn and retraction) is in `JOURNAL.archive.md`, date-ordered, verbatim as of
2026-09-16. Both files live ONLY in this git mirror (since 2026-09-17); edit them here, since
`sync_to_github.sh` no longer copies them. This file keeps only:

0. **Conclusions so far** — the high-level picture, each point backed by the status board.
1. **Status board** — what's been checked and can be trusted, what's retracted, what's open.
2. **Planned analyses.**
3. **Log** — short dated entries (newest first). Put long narrative in the archive, not here.

Labels: **[verified]** = rechecked directly against data, with the check described;
**[reported]** = produced by an earlier session and not independently rechecked since;
**[retracted]** = shown wrong, kept so nobody reuses it.

---

## 0. Conclusions so far (2026-09-17, genome-wide, ~200 donors)

1. **The ONT data is usable at ~30x per haplotype per CpG** once non-CpG junk calls are removed.
   All earlier "sparse coverage" statements were a pipeline bug.
2. **Chemistry (R9.4.1/Guppy vs R10.4.1/Dorado) is real and is not harmonized in the modbeds.**
   - It shifts call-weighted global mCG by ~4 points (R² 0.21–0.26). Measured as an unweighted
     per-CpG mean, the shift nearly vanishes (R² 0.015).
   - The shift is concentrated in low-methylation domains (R10 in-PMD mCG 0.56 vs R9 0.61).
   - It is PC2 of the 10kb matrix (R² 0.30).
   - It does NOT move domain *locations*: R9 vs R10 domain-frequency r = 0.994.
3. **The dominant axis of LCL methylome variation is a continuum of domain hypomethylation, not
   ancestry.**
   - PC1 of the raw 10kb matrix explains 82% of variance. PC1 R²: PMD burden 0.59, HPRC2 QC
     flag 0.22, chemistry 0.19, superpopulation 0.11.
   - In domain bins, donor global state explains R² ≈ 0.91 of between-donor variance, vs 0.34 in
     never-domain bins.
4. **Domain locations are shared across donors ("an extreme version of a common LCL
   methylome").**
   - 41% of 10kb bins are constitutive PMD bins (≥90% of donors); donor-vs-donor PMD-bin
     Jaccard median 0.87.
   - NA19338's PMD bins: 65% constitutive, 15% variable, 15% rare. So it has the common
     domains, deeper, plus a set of rarer ones.
   - The most-methylated donors still show the same domains, just shallowly.
     - HG04187: 10kb mCG 0.685 in constitutive bins vs ~0.83 in never-PMD bins.
     - Its own genome-wide calls have contrast 0.06–0.08.
     - HG02129's own calls carry no contrast at all (~0), although its constitutive bins are
       also lower.
5. **Genome-wide `dnmtools pmd` calls are relative, not a presence/absence phenotype.**
   - Every haplotype gets 1.1–1.9 Gb (median 1.47 Gb; 82% of PMD bp in >1 Mb domains).
   - The inside-vs-outside contrast ranges from −0.005 to 0.36 and tracks global state.
   - Use continuous measures (domain contrast, per-bin mCG) as phenotypes.
6. **Haplotypes agree at domains.**
   - hap1/hap2 PMD-bin Jaccard median 0.91.
   - Mean |hap1−hap2| in 10kb bins: 0.026 in constitutive domains vs 0.016 in never-domain bins.
   - Domain state is shared by both alleles, consistent with the user's WGBS result, so it is
     not a source of ASM.
7. **LCL vs fibroblast PMDs overlap only partly, and the overlap is carried by LARGE domains.**
   - Mean fraction of each fibroblast PMD's bp inside a donor's LCL consensus PMDs, by
     fibroblast PMD size (QC12):

     | size | <100kb | 100–300kb | 300kb–1Mb | 1–3Mb | >3Mb |
     |---|---|---|---|---|---|
     | covered | 0.17 | 0.26 | 0.58 | 0.90 | 0.96 |

   - Fibroblast PMDs over 1 Mb are PMD in ≥95% of donors.
   - Genome-wide hap-consensus Jaccard median 0.51 (all sizes) vs 0.52 (both sets ≥300kb). The
     Jaccard is held down by LCL PMDs covering more sequence (LCL-constitutive domains: 1.08 Gb,
     median 1.2 Mb; fibroblast: 0.94 Gb, median 304 kb), not by the small fibroblast PMDs.
   - Both sets ≥300kb, per donor, medians (hg38 autosomes 2.875 Gb):

     | LCL | fibroblast | intersection | union | Jaccard |
     |---|---|---|---|---|
     | 1.23 Gb (42.9%) | 0.84 Gb (29.3%) | 717 Mb (24.9%) | 1.36 Gb (47.3%) | 0.52 |

     85% of large-fibroblast-PMD bp is in the donor's large LCL PMDs, while only 58% of LCL bp
     is in fibroblast PMDs. LCL-constitutive ≥300kb vs fibroblast ≥300kb: Jaccard 0.55,
     intersection 676 Mb.
   - 65% of LCL-constitutive bins are fibroblast PMDs; 3% of never-domain bins are.
   - Boundaries shared with fibroblasts are gene-enriched (1.13x, permutation p = 0.01);
     LCL-specific recurrent boundaries are not (1.02x, p = 0.17).
8. **Most between-line variance is in domains, and outside regions move with it linearly**
   (QC14, 10kb bins, 201 donors).
   - Share of between-donor variance:

     | bins | share of bins | share of variance | explained by donor state |
     |---|---|---|---|
     | constitutive | 41% | 60% | 76% |
     | rare + variable + common | 28% | 30% | — |
     | never-PMD | 31% | 9.7% | 47% |

   - Donor never-bin mCG = 0.63 + 0.27 × constitutive-bin mCG (chemistry-adjusted, R² 0.84).
     No curvature (p = 0.93); the same slope in the deep, middle and shallow ranges
     (0.30 / 0.26 / 0.33). Outside-domain mCG moves continuously at ~¼ the in-domain rate; there
     is no threshold where it becomes "normal".
   - NB: the per-haplotype "outside own PMDs" value (e.g. NA20762 0.60) is lower than
     never-bin mCG (0.67), because it includes rare/variable bins that are also hypomethylated.
   - Outside domains, state-linked variance sits in intermediate-mCG intergenic bins: never-bin
     R²_state is 0.68 for intergenic vs 0.22 for gene bodies (median mCG 0.88), and 0.67–0.72 for
     bins at 0.6–0.8 mCG vs 0.20 above 0.8. It rises with distance from constitutive domains
     (0.20 at 10–50kb → 0.40 beyond 5 Mb).
   - Open: after regressing out domain depth and chemistry, a residual PC1 carries 27.5% of the
     never-bin variance, and no metadata explains it (superpopulation, sex, passage, era, depth,
     N50: all R² < 0.05).
9. **Metadata (QC13).**
   - Passage: no association (144/155 recorded lines are p5; 46 missing).
   - Age: not available (Coriell lists "Age: No Data").
   - Sex: no domain association. It does associate with the unweighted native global mean
     (R² 0.27), which includes chrX, so that metric shouldn't be used as a covariate.
   - NA- vs HG-prefixed lines (banking-era proxy): slightly deeper domains in NA lines (partial
     R² ≈ 0.02 after chemistry, p ≈ 0.03), but NA lines are mostly R1041 (30/45).
   - Established at Coriell vs externally: nothing.
   - HPRC2 QC flag: strong, because it marks the 3 deepest donors.
   - "Expansion" metrics all collapse onto one axis (Spearman 0.97–0.99): depth, fixed-threshold
     breadth (<0.5, <0.6), relative breadth, and spread into non-constitutive bins. Domains
     deepen in place and neighbouring variable/rare bins deepen with them; there is no separate
     boundary-expansion axis. The metagene boundary-spread metric is pending A01f.
10. **Genetic regulation is not excluded from domains.**
   - HPRC2 promoter mQTLs, among TSS-containing 10kb bins, adjusted for CpG count and TSS count:
     rare/variable-domain bins OR 1.47/1.36, constitutive OR 1.18 (all p < 1e-4, ref = never).
   - Domain-bin variance is dominated by donor state (R² 0.91), but within-donor hap
     differences are larger in absolute terms there (0.00056 vs 0.00021).
   - Genetic signal in domains has to be read against a large state background. The ASM calls
     (C02a) will quantify this directly.

Caveats:
- "Global state" is computed from the same bins, so R²_state is partly circular.
- mQTL tests only cover promoter bins, and S11 lists significant hits only.
- All superpopulation contrasts need chemistry adjustment.

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
- **[retracted] P03 matrix as ASM input.** `results/all_donors/per_sample_chrom/` (4,422 files,
  166G) came from `H02.parse_locus`, which:
  - keeps non-CpG junk calls (HG00097 chr20: 1.85M rows per hap vs 730K CpGs);
  - doesn't merge strands (a minus-strand call sits on the CpG's G);
  - has no flipped-block correction;
  - drops read identity, which the read-level test needs.
  Superseded by C02a; the user can delete the P03 outputs.
- **[verified] New ASM pipeline** (`scripts/asm/`):
  - C01a: shared hg38 region table. 2,019,218 CpG clusters (gap ≤500bp, cut at 1kb, ≥5 CpGs;
    26.1M of 27.7M autosomal CpGs) plus 229 Zink 2018 imprinted DMRs.
  - C02a, per sample × chromosome:
    - same read filters as H02 (≥1 het site in the read span, no HMMFlagger overlap);
    - CpG-checked, strand-merged, flip-corrected calls mapped to hg38;
    - per-read region fractions, then `ont_asm_caller.test_region_reads` (Welch/Mann-Whitney,
      conservative of the two), requiring ≥3 calls per read and ≥3 reads per hap;
    - also writes het-filtered per-CpG counts.
  - C03a: per-chromosome donor-hap × CpG union matrices in hg38 (het-filtered and all-reads).
  - C04a: genome-wide merge. Per-donor genome-wide BH; ASM = q < 0.05 and |Δ| ≥ 0.2; region
    penetrance by chemistry; Zink control by chemistry.
- **[verified] C02a test, HG00097 chr20:**
  - 2.3 min, 2.8 GB; 94% of reads pass the het filter.
  - 49,842 regions tested (median 25 reads per hap); 270 at FDR < 0.05, 192 of them with
    |Δ| ≥ 0.2.
  - 11/12 chr20 Zink DMRs called ASM (GNAS cluster |Δ| 0.67–0.91).
  - CpG output: 730K CpGs per hap, 98.7% on an hg38 CG.
  - `ont_asm_caller` needs `math.comb` (Python 3.8+); C02a shims it for the Python 3.7 allcools
    env.
- **[reported]** Dipcall-style VCFs (G01→G03, jobs 14766207/10/13): validated on 5 donors
  (ts/tv 1.94); full-run completion not rechecked.
- **[reported]** QC09 chain-gap asymmetry: mean 0.6% of the genome (max 1.43%).
- **[reported]** Ancestry PCA: PC1 51%, PC2 18.8%. Output is in `reference/1000G/pca/asm_lr_hprc2/`.

### PMD definition — what "constitutive" means, and alternatives [verified 2026-09-17]

- Current definition: a 10kb bin is **constitutive** if ≥90% of the ~200 donors have it inside
  their own genome-wide `dnmtools pmd` calls (`freq_pmd ≥ 0.9`, B01b). 107,505 bins = 1.08 Gb.
- Alternatives compared (Jaccard vs the current set; donor-level domain depth from each):

  | definition | Gb | Jaccard | corr of donor depth with current |
  |---|---|---|---|
  | freq ≥ 0.9 (current) | 1.08 | 1.00 | 1.00 |
  | freq ≥ 0.8 | 1.15 | 0.94 | — |
  | PMD in all 5 deepest donors | 1.35 | 0.77 | 0.999 |
  | NA20762 (deepest donor) alone | 1.66 | 0.64 | — |
  | NA19338 alone | 1.64 | 0.65 | — |
  | caller-free: mean mCG < 0.70 | 1.76 | 0.59 | 0.990 |

- **Every definition gives essentially the same per-donor phenotype (r ≥ 0.99)**, so conclusions
  about the continuum don't depend on this choice.
- Using the lowest-methylation donors' own PMDs is reasonable (the caller is most reliable where
  contrast is high) but their sets are ~50% larger and include donor-private low regions; the two
  deepest are also the HPRC2-QC-flagged donors. Recommended: keep the frequency-based consensus
  as primary, and use "PMD in all 5 deepest donors" as a maximal-extent sensitivity set.
- A caller-free threshold on mean mCG gives the largest apparent contrast (0.26 vs 0.15) but is
  partly circular (bins are selected for being low in the same data); only use it defined on
  held-out donors.

### Genetic variation and PMDs [in progress 2026-09-17]

- **HPRC2 promoter mQTLs are not depleted in domains** (done, section 0 item 10): among
  TSS-containing bins, adjusted for CpG and TSS count, rare/variable-domain bins OR 1.47/1.36 and
  constitutive OR 1.18 vs never-PMD (p < 1e-4). Limited to promoters and to published
  significant hits.
- **Variant density per 10kb bin** (`B04a_variant_density.py`, job 14778700): SNV / indel /
  ≥50bp-SV counts per donor per bin from G01's per-chrom assembly-vs-hg38 diploid VCFs.
  - chr21 pilot (201 donors): constitutive bins have ~17% more SNVs per donor than never-PMD
    bins (median 16.1 vs 13.7 per 10kb), +3.5 after adjusting for CpG count and gene content;
    indels +0.7; SV +0.10 per bin (median 0 in both).
  - This is germline divergence from hg38, i.e. the known late-replication mutation-rate effect,
    NOT somatic instability in these lines.
- **[verified] The merged per-donor VCFs in `data/vcf/per_donor/` are chr21-only** — G02 ran
  while the G01 array was still going and it tolerates partial chromosome sets. The full
  per-chrom output does exist (`data/vcf_tmp/tmp_<chrom>/<sample>/hap_vs_hg38/diploid.vcf.gz`,
  199–201 donors per chrom). G02/G03 need rerunning (to new filenames) before any cohort VCF is
  used. `data/vcf_tmp` is 1.9 TB, mostly per-chrom `hap1.fa`/`hap2.fa`/`chr*.fa` copies and
  uncompressed `diploid.vcf` — a cleanup candidate once G02/G03 are rerun.
- **Literature expectation** (search 2026-09-17): PMD hypomethylation tracks cumulative cell
  divisions — methylation loss at late-replicating, lamina-associated domains scales with
  population doublings and stops when replication is blocked (Zhou 2018 Nat Genet; Endicott 2022
  Nat Commun). In cancer, PMD depth correlates with somatic mutation density. So the expected
  relationship is: replication timing / mitotic history drives both hypomethylation and elevated
  mutation density, rather than PMDs causing instability. LCL-specific PMD features are
  catalogued in Salhab 2018 (195 methylomes).

### Mechanism, clock, instability, clonality (2026-09-17 evening, all five follow-ups launched)

- **[verified] Replication timing explains the domain map.** `B06a_rt_lad_annotation.py` lifts
  ENCODE/UW Repli-seq wavelet signal for **GM12878** (an LCL — matched cell type) and the UCSC
  laminB1 LAD track (Guelen, fibroblast — the standard cLAD set, not LCL) from hg19 to hg38 10kb
  bins (263,774 bins).
  - Spearman(domain frequency, RT) = **−0.824**; Spearman(mean mCG, RT) = 0.515.
  - Latest-replicating decile: 99.2% constitutive-domain bins, 62% LAD, mean mCG 0.546.
    Earliest decile: 0.3% constitutive, 5% LAD, mean mCG 0.742.
  - LAD fraction by class: never 0.13 → constitutive 0.61.
  - This is the mechanism the literature predicts, measured in our own data.
- **[verified] Solo-WCGW sharpens the phenotype** (`B05a_solo_wcgw.py`, Zhou 2018 mitotic clock).
  3,386,109 solo-WCGW sites in hg38 autosomes (WCGW CpG with no other CpG within 35bp).
  NA19338 hap1: constitutive-domain mCG **15.0%** at solo-WCGW vs 36.6% over all CpGs
  (never-PMD 59.5% vs 67.9%), i.e. depth 44.5 vs 31.2 points. Job 14779380 runs all haplotypes.
- **[verified, chr21 pilot] Variant density is higher in domains** (`B04a`, job 14778700):
  constitutive bins +17% SNVs per donor vs never-PMD bins (16.1 vs 13.7 per 10kb), +3.5 adjusted
  for CpG/gene content; indels +0.7; SV +0.10. QC15 re-tests this genome-wide **with RT as a
  covariate**, which is the question that matters (is the excess just late replication?).
- **[verified] XIST promoter skew now gives real signal** (`scripts/qsub/M03_xist_promoter_skew.py`,
  job 14779488). QC05's whole-gene window gave skew ~0.037; the promoter CpG island gives
  HG00097 0.73/0.29 (skew 0.43) and NA19338 0.57/0.21 (skew 0.36).
  - The window was placed empirically (CpG density peaks TSS−1500→TSS, ~49 CpGs/1.5kb): the local
    `cpgIslandExt.hg38.bed` annotates no island at XIST and its chrX entries look unreliable
    there, while gencode v43 puts XIST at chrX:73,820,649-73,852,723 (−).
- **[verified] RNA marker panel for 200 donors** (`scripts/rna/R01a_rna_markers.py`): HPRC2 Kinnex
  `expression.{plus,minus}.hg38.bw` read **remotely** (no download), 26 genes covering
  proliferation, plasmablast, naive/memory B, activation, senescence and housekeeping, normalised
  to ppm of each donor's total stranded signal.
  - Gotcha: `pyBigWig.stats` must be called with `exact=True`; the default zoom-level
    approximation was off by >10⁶-fold on longer genes (ACTB: −0.33 vs −149,115,297).
- **Genome-wide per-donor VCFs rebuilding**: G02/G03 now write to `data/vcf/per_donor_gw/` and
  `data/vcf/cohort_gw/` (jobs 14779451/14779452), leaving the stale chr21-only outputs untouched.
- **QC15** (`notebooks/qc/QC15_mechanism_clock_instability.ipynb`, job 14779515, held) ties these
  together: RT/LAD vs domains, variant density with RT adjustment, solo-WCGW vs all-CpG depth and
  its chemistry sensitivity, and the QC14 residual axis vs XIST skew and RNA state.
- **C05a** (`scripts/asm/C05a_asm_by_domain.py`, job 14779512, held on C04a) asks whether ASM is
  concentrated in domains and whether it looks genetic (recurrent across donors, mQTL-linked) or
  stochastic (donor-private), plus whether per-donor ASM yield tracks domain depth.

### Jobs (as of 2026-09-17 ~10:30)

Done:
- A01a methcounts: 404/404.
- A01c genome-wide PMDs: 404/404.
- B01a bins: 401/404. HG00272 has no chain; HG01496 hap2 was node-killed and resubmitted as
  14777387.
- B01b/B02a summary (14774195).
- P03 (superseded).

The first A02a bigWig run (14772524) hit the 1h limit in liftOver on every task (0 useful
outputs apart from 6 early files). It was rewritten to use the vectorized chain mapping.

| job | what | depends on |
|---|---|---|
| 14777384 | A02a per-hap hg38 bigWigs `*.hg38.{meth,depth}.bw` (5.5 min, 9.8 GB each) | — |
| 14777387 | B01a rerun, HG01496 hap2 | — |
| 14777427 | C02a ASM calls, 4,444 sample×chrom tasks (chr1 ≈ 8 min) | — |
| 14777428 | C03a donor×CpG union matrices per chrom (hetfilt + all) | C02a, A02a |
| 14777429 | C04a genome-wide ASM merge + Zink/chemistry control | C02a |

**For the user to delete (none of this is used any more):**
```
# old unfiltered methcounts (154G)
rm /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*_hap?.methcounts.tsv.gz
# temp files left by the timed-out liftOver bigWig run (699G)
rm /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*.native.bedGraph \
   /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*.hg38.bedGraph \
   /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*.hg38.sorted.bedGraph \
   /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*.unmapped \
   /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*.chain
# 6 liftOver-era bigWigs (1bp offset on flipped blocks; superseded by *.hg38.meth.bw)
rm /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/data/pmds/*/*_hap?.hg38.bw
# superseded P03 matrix (166G)
rm -r /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/results/all_donors/per_sample_chrom
```

### PMD metagene + size analyses (2026-09-17)

- A01e (`scripts/call_pmds/all_donors/`) builds hg38 reference sets in
  `results/pmd_metagene/ref_sets/`, each also as ≥300kb:
  - lcl_constitutive: 671 domains, 1.08 Gb;
  - lcl_common_plus;
  - fibroblast: 1,410;
  - NA19338 / HG04187 / NA20762 consensus.
- A01f computes metagene profiles from the per-CpG bigWigs with `pyBigWig.stats(nBins)`: 40
  scaled body bins, 200kb flanks in 20 bins per side, flank bins overlapping a neighbouring PMD
  masked. Same design as the igvf YR2 03a–03c notebooks, without window BEDs or tabix.
  About 1 min per hap. Job 14777714, held on A02a.
  - HG00097 hap1 test: flank ~71% vs body ~53% on lcl_constitutive_ge300kb; body depth
    (flank − body) 17.4 there vs 13.1 on fibroblast PMDs.
- `notebooks/qc/QC12_pmd_size_overlap_and_metagene.ipynb`: size-stratified overlap, continuum
  scatter, metagene plots by global mCG and by chemistry. Figures in
  `results/qc/figures/qc12_*`, tables in `results/qc/data/qc12/`. It is re-executed by
  U01c once all profiles exist.

### UCSC track hub [verified built]

- **Gradient hub:** `/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2_gradient`, built by
  `scripts/ucsc/U01b_build_gradient_hub.py`.
  - 20 R941 donors' hap1 per-CpG mCG, evenly spaced in call-weighted global mCG (NA19338 0.531
    → HG02129 0.748), ordered and viridis-coloured. Selection in `hg38/selected_donors.tsv`.
  - R941 only, to keep chemistry out of the gradient. The low end is sparse: NA19338 is the only
    R941 donor below 0.60.
  - Reference tracks: fibroblast PMDs, LCL constitutive domains, PMD frequency, mean mCG.
  - 9 donors' bigWigs weren't written yet at first build. Job 14778024 (U01c, held on A02a +
    A01f) rebuilds both hubs and re-executes QC12.
- **Gradient hub, revised (2026-09-17 evening):**
  - Now both chemistries, per user request, to smooth the low end. The two deep outliers
    (NA20762 0.415, NA19338 0.531) are kept, and the other 18 picks are evenly spaced from
    0.567 to 0.748 (HG02129). Chemistry is in each track label; files are named by sample.
  - HG02165 (rank 2 in the first, R941-only build) is 15th lowest globally out of 201
    (9th lowest in-PMD). HG02129 is highest globally and in constitutive bins (201/201).
  - Unreferenced files from earlier selections can be deleted:
    `cd /u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2_gradient/hg38 && rm [0-9][0-9]_*_hap1_mCG.bw HG00290_hap1_mCG.bw HG00642_hap1_mCG.bw HG00658_hap1_mCG.bw HG01928_hap1_mCG.bw HG02004_hap1_mCG.bw HG02273_hap1_mCG.bw HG02809_hap1_mCG.bw HG03583_hap1_mCG.bw NA18952_hap1_mCG.bw NA18976_hap1_mCG.bw NA19036_hap1_mCG.bw NA20346_hap1_mCG.bw NA21093_hap1_mCG.bw NA21144_hap1_mCG.bw`
    (check `trackDb.txt` first if the selection has changed again).
- **Why the chr1 centromere / 1q12 is empty in the browser [verified]:** it is an hg38
  limitation, not coverage.
  - hg38 chr1:126–142 Mb is 100% N; 121–125 Mb is modeled satellite with CpGs but no chain
    alignment.
  - In the assemblies the region is present: HG00097 hap1 chr1 is one 251.6 Mb contig with
    24.9 Mb between the last p-arm and first q-arm hg38 blocks. That stretch has 26–98K CpGs per
    2 Mb and ~900–1,050 read starts per 2 Mb (genome average ~1,200), median read span ~37 kb.
  - NA19338 hap1 chr1 is split into two contigs at the centromere, with ~18 Mb of unaligned
    satellite at the p-arm contig end, also covered.
  - This unmapped satellite is a large part of the 18% of CpGs lost in hg38 mapping. A CHM13
    view (`*_vs_CHM13.chain.gz` exist) would show these regions.
- **Main hub:**

- Location: `/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2`, same layout as the lab's other hubs.
  Built by `scripts/ucsc/U01a_build_track_hub.py` (re-runnable), 967 MB.
- superTrack `LCL_population`: fibroblast PMDs, and 10kb mean mCG / SD / PMD frequency /
  caller-free low-mCG frequency / PMD-boundary frequency.
- superTracks `NA19338` and `HG04187` (lowest and highest global mCG): hap1, hap2 and consensus
  PMDs, plus per-CpG mCG per hap. The per-CpG tracks are added when the script is rerun after
  A02a finishes.
- superTrack `LCL_donors_10kb`: all 401 haps at 10kb (hidden by default; blue R941, orange R1041).
- Per-CpG bigWigs for everyone (~290 GB) were deliberately not copied.

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

**A–D are run genome-wide** (results summarized in section 0). Their outputs:
- `results/meth_bins/` (B01b)
- `results/meth_bins/boundaries/` (B02a)
- `results/meth_bins/qc_genetics/` (B03a: PMD QC, hap agreement, mQTL-by-domain logit,
  per-bin variance decomposition)

**A–C scripts** in `scripts/meth_bins/`:
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

**Next steps for the PMD arm (recommended order, 2026-09-17):**
1. **Solo-WCGW mitotic clock.** Zhou 2018's sharpest PMD phenotype is mCG at solo-WCGW CpGs
   (a W-C-G-W CpG with no other CpG within 35bp). Cheap here: annotate them from hg38, recompute
   per-donor domain depth on that subset, and check whether it tightens the continuum and its
   correlation with anything (it is the standard proxy for cumulative divisions).
2. **Replication timing / LAD overlap.** Public LCL (GM12878) Repli-seq and LADs are the direct
   test of the mechanism, and would explain both the domain map and the SNV-density result.
3. **Finish the variant work**: genome-wide B04a, then rerun G02/G03 for genome-wide per-donor
   VCFs, and ask whether domain depth associates with a donor's own variant burden.
4. **ASM inside domains** (needs C02a/C04a): domains are where cross-donor variance lives, so
   the interesting question is whether allelic asymmetry is also concentrated there, and whether
   ASM in domains is genetic (mQTL-like) or stochastic/epigenetic drift.
5. **Explain the residual axis** from QC14 (27.5% of outside-domain variance, no metadata
   correlate): candidates are LCL clonality (the XIST/XCI metric once fixed to the promoter),
   EBV copy number, and B-cell differentiation state from the Kinnex RNA.

**Also pending:**
- Genome-wide version of the chr20 boundary/gene analysis (A03d) once 14772523 lands.
- Recompute QC10.
- Chemistry-stratified replication of every per-donor result above.

---

## 3. Log (newest first)

### 2026-09-17 (late evening) — all five PMD follow-ups launched
- RT/LAD annotation built; domain frequency vs replication timing Spearman −0.82 (mechanism
  confirmed in our data).
- Solo-WCGW sites annotated (3.39M) and per-haplotype job running; pilot shows much deeper
  domain hypomethylation (15% vs 37%).
- Variant-density array running; chr21 pilot shows +17% SNVs in domains.
- XIST promoter skew fixed and running; RNA marker panel done for 200 donors (remote bigWig).
- G02/G03 rebuilding genome-wide per-donor VCFs.
- QC15 and C05a written and held on their inputs.

### 2026-09-17 (night)
- ≥300kb LCL-vs-fibroblast Jaccard + genome coverage.
- Gradient hub rebuilt with both chemistries and outlier-aware spacing.
- chr1 centromere gap diagnosed (hg38, not coverage).
- QC13 (`notebooks/qc/QC13_pmd_expansion_metadata.ipynb`): PMD expansion metrics vs metadata.
  Nothing but chemistry and the HPRC2 QC flag; a weak NA-prefix signal; no age data.
- QC14 (`notebooks/qc/QC14_variance_inside_outside_pmds.ipynb`, genome-wide successor of the
  chr20 QC11): variance share by domain class, outside-vs-inside coupling (linear, slope 0.27),
  and where outside-domain state variance lives. Also found an unexplained residual axis.
- ASM calls 1,401/4,444; bigWigs 402/404 (HG00272 has no chain).

### 2026-09-17 (evening)
- User hypothesis confirmed: big fibroblast PMDs are conserved in LCLs (>1 Mb: 90–96% covered),
  small ones mostly aren't. Jaccard is insensitive to the size cut.
- HG04187 has shallow but present domains.
- Added PMD reference sets (A01e), bigWig metagene (A01f, job 14777714), QC12 notebook, the
  gradient hub (U01b), and a held rebuild job (U01c, 14778024).

### 2026-09-17 (later)
- Genome-wide B01b/B02a results in.
- Found and fixed:
  - the A02a liftOver timeouts (vectorized mapping, new `*.hg38.{meth,depth}.bw`);
  - the P03 ASM-input problems (junk calls, unmerged strands, no flip correction, no read
    identity), replaced by `scripts/asm/` C01a–C04a (tested on HG00097 chr20) and submitted.
- The QC04 vs B01a global-mCG chemistry discrepancy is a metric difference (call-weighted vs
  per-CpG mean), not an error.
- B03a: genome-wide PMD QC, hap agreement, mQTL-by-domain, variance decomposition.
- Built the UCSC hub. Moved the journal to the mirror only. Added section 0 (conclusions).
- Pending: Zink imprinted-DMR ASM by chemistry (C04a), bigWigs, CpG matrices, hub per-CpG
  tracks.

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
