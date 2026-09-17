# Journal

Reverse-chronological. One entry per session that changes project state. This file — not a
proliferating set of dated `md/YYYYMMDD_*.md` files — is the single place a new session should
read first. See `README.md` for current project state/design; this file is *why* it got that way
and what's still open, in the order it happened.

Entry template:

```
## YYYY-MM-DD — one-line summary

**State at start:** …
**Decisions made:** …
**Produced:** (files/commits)
**Open / next:** …
**If resuming, read:** (the one or two files that matter most)
```

---

## 2026-09-16 — data/ reorg broke every hardcoded PROJDIR-relative path; found via H01/H02 audit, fixed, HMMFlagger gap discovered and backfilled, pipeline resubmitted end-to-end

**State at start:** job 14756945 (genome-wide H01, submitted 2026-09-15) had finished. User had
manually moved its output (`scripts/harmonize_hg38/pilot/tmp_<chrom>/<sample>/het_in_hap{1,2}.bed`)
to a new home, `data/het_snps/`, as part of a broader reorg that (unnoticed until this session)
had already relocated `modbed/`, `chains/`, `assemblies/`, `hmmflagger/` from `PROJDIR` root to
`PROJDIR/data/` sometime this week. User asked whether `P03_build_donor_chrom_matrix.py` /
`P03_build_donor_chrom_matrix_array.sh`'s paths were still correct.

**Decisions made:**
- **They weren't — and not just the het-sites path.** `P03_build_donor_chrom_matrix.py` imports
  `load_het_sites`/`load_hmmflagger`/`parse_locus` from `H02_filtered_locus_matrix.py` rather than
  defining its own paths, and that file had **four** stale constants left over from before the
  `data/` reorg: `MODBED_DIR`, `CHAIN_DIR`, `HMM_DIR` (all `PROJDIR/<name>`, should be
  `PROJDIR/data/<name>`) and `load_het_sites()`'s hardcoded `HERE/tmp_<chrom>/<sample>/...` (`HERE`
  = the pilot script dir — the pre-move location, not `data/het_snps/` where the files actually
  are now). Same four-way staleness also present in `H01_call_hap_het_sites.py` (`ASSEMBLY_DIR`,
  `CHAIN_DIR`, plus its own `HERE`-based output dir), `G01_call_hap_vs_hg38.py` (`ASSEMBLY_DIR`),
  `P04_read_length_coverage_qc.py` (`MODBED_DIR`, `CHAIN_DIR`), and three download scripts
  (`A01a_download_hprc2_files.sh`: `MODBED_DIR`, `ASSEMBLY_DIR`; `A01c_download_liftover_chains.sh`
  + its `test/` variant: `CHAIN_DIR`). Fixed all of it to point at `data/`; `H01`'s output (and
  `H01_call_hap_het_sites_array.sh`'s matching skip-check) now goes straight to
  `data/het_snps/tmp_<chrom>/<sample>/`, so a future rerun never needs another manual move.
- **This wasn't hypothetical — it already caused real, silent-ish data loss in job 14756945.**
  Checked the 4444-task output: only 4364/4444 (sample, chrom) pairs actually produced
  `het_in_hap1.bed` (80 empty output dirs). Root-caused via the task logs: **HG02280 and HG00272
  failed all 22 of their chromosomes** — `samtools faidx`/chain-file `FileNotFoundError` on the
  exact stale `PROJDIR/assemblies` and `PROJDIR/chains` paths just fixed (44 of the 80 pairs,
  explained). The remaining ~36 scattered pairs (weighted toward the largest chroms — chr5:8,
  chr3:8, chr1:7) show no traceback, just truncated minimap2 output — consistent with the
  `h_rt=1:00:00` timing margin the 2026-09-15 entry already flagged as unverified at chr1-3 scale,
  not the path bug. **Resubmitted the full H01 array** (job **14764732**) rather than hand-picking
  80 task IDs — skip-if-exists means the 4364 good pairs no-op instantly, only the 80 gaps
  actually rerun, now against the corrected paths.
- **HMMFlagger track: only 1/202 donors had it** (`HG00097`, fetched by hand during pilot
  development) — `data/hmmflagger/` was effectively empty at the donor-panel scale. Confirmed
  this is per-donor, not per-chrom (`HMMFlagger.ONT.bed.gz`, one file per sample, same
  `hprc-epigenome/samples/<SAMPLE>/` layout `A01a` already uses for modbed) — so **202 files
  needed, not 4444**, and unlike `A01a`'s multi-GB modbed/assembly downloads these are tiny
  (~3.5KB observed for HG00099). `A01a` never fetched it (it only pulls modbed + assembly) — this
  was a genuine gap, not a move-related regression. Wrote `A01d_download_hmmflagger.sh` (+
  `test/` variant, matching `A01a`/`A01c` conventions exactly: manifest-row-indexed SGE array,
  skip-if-exists, `.partial`-then-`mv`) and submitted it, **job 14764729**.
- Why this matters more than an ordinary missing-file bug: `load_hmmflagger(sample)` runs
  **unconditionally**, before P03's per-hap `try/except FileNotFoundError` (that block only
  wraps `load_het_sites`). So every one of the 4444 P03 tasks would have hit a hard
  `FileNotFoundError` and failed outright rather than quietly under-counting — loud, not silent,
  but still would have burned the entire `-tc 60` throttle failing instantly if submitted before
  this was caught.
- **Submitted the actual production run**, chained rather than run by hand in sequence: added
  `A01d_download_hmmflagger` to `P03_build_donor_chrom_matrix_array.sh`'s `-hold_jid` (now
  `H01_het_sites_array,A01d_download_hmmflagger`) so P03 can't start on a donor before both its
  het-sites and its HMMFlagger track exist, then submitted P03 itself — **job 14764741**, held
  until 14764732 and 14764729 both finish.

**Produced:** path fixes in `H01_call_hap_het_sites.py`, `H02_filtered_locus_matrix.py`,
`G01_call_hap_vs_hg38.py`, `P04_read_length_coverage_qc.py`,
`scripts/download/A01a_download_hprc2_files.sh`, `scripts/download/A01c_download_liftover_chains.sh`
(+ `test/` variant), `scripts/harmonize_hg38/all_donors/H01_call_hap_het_sites_array.sh`;
new `scripts/download/A01d_download_hmmflagger.sh` (+ `test/` variant); updated `-hold_jid` in
`P03_build_donor_chrom_matrix_array.sh`. Jobs: 14764732 (H01 gap-fill), 14764729 (HMMFlagger
backfill), 14764741 (production P03, held on both) — all in progress as of this entry.

**Open / next:**
1. **Check all three jobs next session** — `qstat -u terencew`, then
   `grep -a real logs/P03_donor_chrom_matrix_array.*` once P03 actually starts running (it's
   sitting in `hqw` until the two holds clear) to get the first real chr1-scale timing number the
   2026-09-15 entry flagged as never measured — raise `h_rt` before the full 4444 run if a chr1
   task is close to the current 1-hour limit, since that's the same margin that plausibly caused
   this session's ~36 scattered H01 timeouts.
2. Once P03 finishes, spot-check `results/all_donors/per_sample_chrom/` coverage against the
   4444-row task list the same way this session did for H01 — confirm no analogous silent gap
   before treating the matrix as complete.
3. **The purity/k-filter finding from 2026-09-15 is not a bug and won't be fixed by this
   pipeline running cleanly** — worth restating here since it's easy to lose track of once the
   plumbing issues are resolved: `k=1` het-site filtering removes reads that *can't* carry
   phasing information, but H03 found raising `k` past 1 buys essentially no purity gain (≤0.08
   change across `k=0..30`), and H04 found only 32.7% of loci are even genuinely bimodal in the
   first place — the rest of the residual per-read disagreement needs read
   sequence/genotype information the modbed format doesn't carry, not a better filter. This is a
   real ceiling to disclose in methods/limitations, not something the genome-wide run will
   quietly resolve.
4. `scripts/github/sync_to_github.sh` still doesn't rsync `scripts/harmonize_hg38/` (carried over
   from 2026-09-15, unaddressed this session too) — all the path fixes above live only in the
   working dir until that's fixed and run.

**If resuming, read:** this entry, then check job status first (`qstat -u terencew`), then
`scripts/harmonize_hg38/pilot/H02_filtered_locus_matrix.py` to confirm the path fixes are still
in place before trusting any P03 output.

---

## 2026-09-16 — same day, later: population-genetics QC sweep (coverage, methylation, PCA, SNP/CpG landscape, dipcall VCFs) — found and confirmed a real ancestry-driven detection-rate confound, documented upstream in ont_asm_caller

**State at start:** the entry above had just fixed the post-reorg path bugs and resubmitted
H01/P03. User asked for a broad QC sweep motivated by the README's new framing (population-scale
sequence determinants of ASM, not just detection): coverage QC, global methylation vs.
covariates, LCL monoclonality literature, XCI feasibility, SNP/CpG landscape, PMD feasibility,
dipcall-style per-donor VCFs, ancestry PCA, and population-genetics stats (AF/LD). Dispatched as
a series of parallel forks plus direct work; results below are what actually finished.

**Decisions made / findings:**
- **Coverage QC** (`notebooks/qc/QC03_actual_vs_reported_ont_coverage.ipynb`): HPRC2's reported
  `coverage_ont` (Supp Table S6) is diploid/total, not per-haplotype — measured actual
  hap1+hap2 depth vs. reported, median ratio 0.94 (expected shortfall: modbed only contains
  mapped+haplotype-assigned reads). **Recommend reporting both diploid and per-haplotype depth**
  in any donor QC table — diploid matches HPRC2's own convention, but per-haplotype is what
  actually gates power in H02's het-site filter. Two donors (HG02293, NA19776) flagged at
  ratio 0.66, well below the rest (0.86–1.18) — only single-window-sampled, needs a genome-wide
  check before treating as real.
- **Global methylation vs. covariates** (`QC04_global_methylation_covariates.ipynb`, full
  219-donor cohort, 438 hap-observations, `scripts/qsub/M01_global_methylation_array.sh` job
  14765944): range 0.521–0.726. Directly quantified: **ancestry alone R²=0.115** (F=13.6,
  p=2×10⁻¹⁰), sex alone 0.017, haplotype (within-donor) 0.002, all three combined only 0.125 —
  ancestry is doing essentially all the explanatory work in the "common covariates" set, and it
  still leaves ~87.5% of variance unexplained. **Passage number** (found in
  `reference/hprc2/hprc2_supp.xlsx` Table S15 — missed by the first covariate-scoping pass,
  which only checked the manifest + public 1000G metadata) looked like a real additional
  covariate at small n (15 donors: R² 0.477→0.647 adding passage, p=0.0056) but **did not hold
  up at full cohort scale** (R² 0.161→0.168, p=0.099) — passage is confounded with ancestry
  cohort-wide (passage=3 is 12/12 EUR), so the small-n result was likely ancestry leaking through
  as apparent passage effect. Age and passage-8-12 donors remain genuinely unavailable (no
  age/DOB field anywhere public for this cohort; passage 8-12 donors lack modbed data entirely).
- **XIST skew QC** (`QC05_xist_skew_qc.ipynb`, 98 female donors, coordinates verified via
  Ensembl REST): skew is small and unimodal (mean 0.037, max 0.117) — **not** the bimodal
  clonality signature a working assay should show. Diagnosed, not treated as a real negative:
  it measured the whole ~32kb XIST gene body instead of the diagnostic 5' promoter/CpG island,
  which would dilute any real skew toward exactly this muted pattern regardless of true
  clonality. Narrowing the window is the next step, not a pipeline rebuild. Also: even a clean
  measurement here would only rule out severe/skewed monoclonality, not a balanced 2-clone
  mixture — XCI skew is a chrX-specific readout, not a substitute for an autosomal clonality
  metric.
- **LCL monoclonality literature** (no code, WebSearch): confirmed via Plagnol et al. 2008
  (*PLoS ONE*, PMC2494943) — using XCI skew across 1,174 LCLs, only 52–68% had balanced XCI vs.
  88–92% in matched peripheral blood, ~60% pauciclonal, **≥22% effectively monoclonal**. Real,
  quantified confound for bulk LCL-derived ASM, not hypothetical. Cheapest workaround: the XIST
  skew metric above, once fixed to the promoter window, doubles as a per-donor clonality QC
  filter/covariate.
- **PMD feasibility**: confirmed in literature that EBV-transformed LCLs specifically (not just
  cancer generally) reliably produce PMDs (Frontiers in Genetics 10.3389/fgene.2017.00076 and
  corroborating sources) — this cohort should be treated as PMD-positive by strong prior. No
  real computational first pass yet; drafted (unexecuted) skeleton at
  `notebooks/qc/QC06_pmd_windowed_methylation_DRAFT.ipynb` (windowed mean methylation +
  CpG-density correction + segmentation).
- **HPRC2 paper mining**: the actual paper (bioRxiv 2026.07.21.739710v1, PMC13419748) and its
  22-table supplementary xlsx were already sitting locally, unopened
  (`reference/hprc2/hprc2_supp.xlsx`) — opened directly rather than re-fetching. Table S11
  ("Significant promoter mQTL with lead variants," 80,854 rows) is HPRC2's own var-CpG table
  (promoter-scoped, not genome-wide) — directly reusable as a positive control. Table S15 is
  where passage number and karyotype live (see above); only 4/234 donors show non-standard
  karyotype, all benign constitutional variants, not culture-induced instability.
- **Dipcall-style per-donor VCFs**: `G01_call_hap_vs_hg38.py` (dipcall-style hap-vs-hg38 calling,
  already existed from 2026-09-15 but had never actually been run and had the same stale-path
  bug as everything else — fixed) → new `G02_build_donor_snp_vcf.py` (biallelic-SNP filtering)
  → new `G03_merge_cohort_vcf.py` (cohort-merged **and** per-superpopulation-merged VCFs).
  Validated on 5 donors (one per superpopulation): 75–88K variants/donor, cohort-merge ts/tv=1.94
  (healthy). Chained full run submitted: G01 array job **14766207** (4444 tasks) →
  G02 job **14766210** → G03 job **14766213**, still running as of this entry. Confirmed asm_lr's
  original dipcall tooling (`D03`/`D04`, real `dipcall-aux.js`) was the precedent, but G01's
  direct-phasing approach is correct for this input shape (see 2026-09-15 entry for why
  `vcfpair` doesn't fit an assembly-vs-reference comparison) — not a regression to fix.
- **Known permanent gap**: HG00272 has no chain file on HPRC2's bucket at all (confirmed 404,
  not a download bug) — that donor's liftover-dependent outputs (H01, G01/G02) will have a
  contained, permanent gap unless HPRC2 deposits it or chains are derived independently.
- **Ancestry PCA**: built a reusable pipeline (now living in the shared
  `/u/project/cluo/terencew/claude/reference/1000G/scripts/run_1000g_pca.sh`, since this is
  useful across projects, not just this one) — subsets the 1000G high-coverage (NYGC, hg38)
  panel to an arbitrary donor list, biallelic SNPs, MAF≥0.05, LD-pruned, `plink --pca`. Run on
  221/229 HPRC2 donors (chr1-3, 104K pruned SNPs): **PC1 51.0%, PC2 18.8%, PC3 8.7%** of
  variance — real, strong continental structure, and zero label mismatches between the 1000G
  pedigree's superpopulation call and this project's own manifest. Output and notebook moved to
  `reference/1000G/pca/asm_lr_hprc2/` and `reference/1000G/ipynb/PCA01_1000g_ancestry_pca.ipynb`
  per the same "shared reference, not project-specific" reasoning.
- **Population genetics stats** (`QC08_popgen_af_ld_hetsite_overlap.ipynb`): using the PCA's
  chr1-3 1000G-subset VCFs. Real ancestry divergence even among panel-common (MAF≥5%) sites —
  EAS: 14.6% become monomorphic and 25.1% become rare within-population vs. AFR's 3.1%/15.4%.
  LD decay is textbook ancestry-structured in this exact panel: AFR r²=0.40 at 0-10kb vs. EAS
  0.78. Attempted het-site (H01) vs. 1000G-panel overlap (>99% no-match) but flagged its own
  caveat: the AF lookup used was already MAF-pruned for the PCA, so "no match" can't yet
  distinguish "genuinely novel" from "present but rare" — needs a rerun against an unfiltered AF
  source. **Decided against re-deriving an unfiltered 1000G table for this — gnomAD is the
  better source** (order-of-magnitude larger per-ancestry N, ships genome-wide per-population
  AF/AC/AN tables as its primary product) for AF-based rare/common classification specifically,
  while 1000G/NYGC remains right for PCA/LD (needs named per-sample genotypes gnomAD doesn't
  publish). Not yet downloaded — plan is to query gnomAD's public GCS-hosted site VCFs by region
  (remote tabix, no bulk download) once there's a concrete het-site position list to look up, not
  pull the whole dataset onto an already ~95%-full filesystem.
- **The headline finding: ancestry-driven heterozygosity density is a real, large, measured
  confound.** Computed directly from H01's own output, 197 donors with complete data:

  | superpop | n donors | mean het sites/Mb |
  |---|---|---|
  | AFR | 54 | 1470.4 |
  | SAS | 36 | 1166.2 |
  | AMR | 40 | 1137.3 |
  | EUR | 29 | 1130.0 |
  | EAS | 38 | 1073.4 |

  One-way ANOVA: **R²=0.917, F=530, p=1.6×10⁻¹⁰²** — superpopulation alone explains 91.7% of
  variance in het-site density. AFR donors carry ~37% more heterozygous sites/Mb than EAS. Since
  every het-site-filtered test (H02's `k≥1`, and any ASM test that needs a phasing-informative
  read) can only fire on reads spanning a het site, this means **raw cross-ancestry "ASM
  detection rate" comparisons are confounded by this mechanical effect before any real biology
  is considered** — a population with lower baseline heterozygosity will show fewer detectable
  loci even if the true underlying regulation is identical.
- **Decided this belongs in the caller, not the filter, and documented it as such in
  `github/ont_asm_caller`** (a separate repo/project, not just this one) — see that repo's
  2026-09-16 JOURNAL/PRIORITIES-item-8/README entries, pushed as commits `3a4d202` and
  `412069e`. Explicitly rejected the tempting fix (ancestry-conditional filter thresholds) as
  worse than the problem — it trades a quantifiable confound for an ancestry-conditional
  analytical choice. The intended fix is a caller-level output change: report per-locus
  het-site-informative-read count as an explicit power/confidence indicator (extending the
  caller's existing measured-nuisance-parameter posture for DE/dispersion/call-error), so
  low-heterozygosity-population loci come out "underpowered" rather than silently "not detected."
- **Housekeeping**: two forks briefly duplicated work on the same methylation-covariates
  question after a mid-session file rename wasn't communicated to an already-running fork —
  caught and cleaned up (one stale `QC03_global_methylation_covariates.ipynb` deleted, canonical
  version is `QC04`). `notebooks/qc/` numbering: QC01-QC08 now assigned, QC06 is an unexecuted
  draft, QC07 (SNP/CpG landscape: het-SNP/CpG counts, CpG-disrupting SNP rate, CpG-density
  architecture, all ancestry-stratified) was still running as of this entry.

**Produced:** `notebooks/qc/QC03`–`QC08` (QC06 draft, QC07 still pending), `scripts/qsub/
M01_global_methylation_array.sh`, `scripts/harmonize_hg38/all_donors/G02_build_donor_snp_vcf.py`,
`G03_merge_cohort_vcf.py` + their array/qsub scripts, `results/qc/data/het_site_density_by_superpop.tsv`,
`reference/1000G/{scripts,pca/asm_lr_hprc2,ipynb,tsv/meta}/*` (shared reference dir, not this
project). Jobs in flight: 14765944 (done, M01), 14766207/10/13 (G01→G02→G03, running).
`github/ont_asm_caller` commits `3a4d202`, `412069e` (pushed).

**Open / next:**
1. **QC07 (SNP/CpG landscape + CpG-architecture) still running** — check on completion; it's the
   direct input to "do we gain/lose CpGs per ancestry" and the corrected CpG-disrupting-SNP rate.
2. **G01→G02→G03 dipcall VCF pipeline still running** (4444-task G01 array) — once done, gives a
   native (non-1000G-borrowed) per-donor and per-superpopulation genotype/AF source.
3. **Rare/private-variant extension planned, not yet built**: per-donor private CpG/SNP fraction
   by superpopulation, and a corrected common-vs-rare het-site classification against gnomAD
   (not the MAF-pruned 1000G table QC08 used). Blocked on QC07/QC08 informing exact scope, and
   on deciding gnomAD access approach (remote region query vs. download).
4. **The actual detection-rate effect of the het-site-density confound on real ASM calls** is
   still unmeasured — needs the beta-binomial caller actually run on this cohort's filtered
   matrix (once P03 finishes) stratified by ancestry, compared against the het-density numbers
   above.
5. Carried over from the entry above: chr1-scale P03 timing came in at 24-35 min (within the
   1hr `h_rt` budget, no change needed) — confirmed this session, no longer open.
6. `scripts/github/sync_to_github.sh` still doesn't cover `scripts/download/` or
   `scripts/harmonize_hg38/all_donors/` (carried over, unaddressed again this session).

**If resuming, read:** this entry, then `github/ont_asm_caller/JOURNAL.md`'s 2026-09-16 entries
for the cross-project confound (it affects both codebases), then check `qstat -u terencew` for
the still-running QC07/G01-G03 work before trusting anything downstream of them.

---

## 2026-09-16 — same day, later still: PMDs confirmed directly in this data, haplotype-asymmetry checks (corrected), read-length/HPRC2-supp/reference-bias/gnomAD questions answered

**State at start:** the entry above had just found and documented the ancestry-driven het-site
density confound. This entry covers the follow-up QC work and conceptual questions from the
rest of the session.

**Decisions made / findings:**
- **ONT read length**: no need to recompute — already in HPRC2's own Supp Table S6
  (`results/qc/data/supp_seq_qc.csv`, `read_N50_ont`). Median N50 81.4kb, range 33.6–115.3kb.
  Unlike het-site density, **read length shows minimal ancestry stratification** (AFR 84.2kb to
  AMR 79.2kb, ~6% spread) — reassuring that this particular technical variable isn't feeding the
  ancestry confound the way variant density is. R1041 chemistry slightly longer than R941
  (83.5kb vs 80.5kb median), already tracked elsewhere in this project's docs.
- **Full HPRC2 supplementary xlsx inventory** (22 tables, `reference/hprc2/hprc2_supp.xlsx`) —
  checked directly rather than assumed. **S10/S11 (the var-CpG/promoter-mQTL tables) have no
  superpopulation or population column at all** — both are pooled across all donors. So a
  superpopulation-stratified var-CpG analysis is NOT recoverable from HPRC2's own published
  tables; it would have to be built from scratch (join S11's lead-variant positions against our
  own per-superpopulation methylation data, or derive fresh from the G01-G03 VCFs).
- **Reference bias / WASP discussion** (conceptual, no new code): this project's per-haplotype
  own-assembly mapping design (hap1 reads → hap1's own assembly, hap2 reads → hap2's own
  assembly) structurally avoids the classical reference-allele-favoring bias WASP corrects for,
  since neither haplotype plays the role of "the reference" the other is penalized against.
  H01/H02's `k≥1` het-site filter is the structural analog of WASP's read-remapping check —
  "only trust reads whose haplotype provenance is actually verifiable" — just implemented via
  differentiating-site presence instead of allele-flip-and-remap. Residual concerns, not
  eliminated by this design: liftover completeness can still differ by haplotype (see chain-gap
  findings below), and HMMFlagger's assembly-reliability exclusion plays an analogous
  "don't trust unreliable regions" role.
- **gnomAD**: confirmed no local copy exists anywhere under `/u/project/cluo` (checked directly).
  Decided against bulk-downloading it — full site VCFs run hundreds of GB to low TB, and the
  filesystem is already ~95% full. Recommended approach for the rare-variant follow-up: query
  gnomAD's public BigQuery dataset or its GCS-hosted site VCFs by region (remote, no local
  storage) once there's a concrete het-site position list to look up — `bq`/`gcloud`/`gsutil` are
  already installed and authenticated in this environment. gnomAD is the right AF source
  specifically for rare/common classification (much larger per-ancestry N than 1000G); 1000G/NYGC
  remains right for PCA/LD (needs named individual-level genotypes gnomAD doesn't publish).
- **PMDs — confirmed directly in this data, not just the literature prior.** Real executed check
  (`notebooks/qc/QC06_pmd_windowed_methylation.ipynb`, replacing the earlier unexecuted draft),
  chr20, 3 donors spanning the global-methylation range, 20kb windows with a CpG-density filter
  (barely binding — median 800-900 calls/window, this cohort's ONT coverage is dense enough) and
  real segmentation (≥5 consecutive low-methylation windows = domain, explicitly excluding
  single-window CpG-island/promoter dips):

  | donor (global meth) | % of chr20 in PMD-like domains | largest domain |
  |---|---|---|
  | NA19338 (0.523, low tail) | **78.3%** | up to 4.36 Mb |
  | HG01981 (0.648, median) | 69.1% | up to 4.72 Mb |
  | HG04187 (0.723, high) | **22.4%** | up to 900 kb |

  The low-global-methylation donor shows both more domain coverage AND larger individual domains
  than the high-methylation donor — the actual PMD signature (large, sustained low/intermediate
  domains), not just a uniformly shifted mean. This meaningfully raises confidence beyond the
  literature prior alone. Caveat: single chromosome, 3 donors, threshold-based (not HMM)
  segmentation — exact percentages are threshold-sensitive, the donor-to-donor gradient is not.
- **Haplotype-asymmetry checks** (`notebooks/qc/QC09_haplotype_asymmetry_chain_gaps.ipynb`, 15
  donors, 3/superpopulation), using existing chain files — no new alignment needed:
  - SV-scale (≥50bp) gaps: mean **5.75% of the genome** sits in gapped/structurally-divergent
    regions per haplotype (range 4.49-6.90%). AFR haplotypes average ~9% more SV-scale gaps than
    EUR — same direction as the het-SNP density finding, much smaller magnitude.
  - Within-donor hap1-vs-hap2 asymmetry in *total gapped bases* averages 10.3% (up to 25.9%,
    HG01891) — but that 10.3%/25.9% is a **relative** difference between two already-modest
    (~5-6%-of-genome) quantities, not a fraction of the genome. **Re-expressed as % of genome,
    the actual hap1-vs-hap2 imbalance is much smaller: mean 0.6%, max 1.43%** (HG00438) and 1.34%
    (HG01891) — caught and corrected this session after the first summary conflated the two
    denominators; worth remembering when this number gets cited again.
  - Post-liftover CpG-count balance (8 donors with available P03 output — all EUR/EAS/AMR, P03
    hadn't reached AFR/SAS yet): tiny, 0.07-0.28% hap1-vs-hap2 difference, essentially at parity.
  - Cross-check: correlation between a donor's chain-gap asymmetry and its post-liftover CpG
    imbalance is **-0.36 (n=8)** — weak, wrong-signed, no support for a detectable downstream
    effect. **Incomplete, not conclusive**: this is only tested in EUR/EAS/AMR, exactly where the
    asymmetry is smallest — needs rerunning once P03 reaches AFR/SAS (where HG01891's 25.9%/1.34%
    outlier lives).
- **Dipcall pipeline (G01→G02→G03) confirmed working**, no further findings beyond the earlier
  validation (5 donors, ts/tv=1.94) — full-scale run (jobs 14766207/10/13) still in progress.
- **QC07 (SNP/CpG landscape + CpG architecture) still running as of this entry** — genuinely
  long-running (multiple hours), no results yet to report.
- **User request for next session**: run an HMM-based PMD caller (upgrade from this session's
  threshold/windowed-segmentation first pass in QC06 — a real HMM, e.g. MethylSeekR-style or
  methylpy's PMD caller, would replace the "≥5 consecutive low-methylation windows" heuristic
  with a proper state-based segmentation) and generate BigWig files from these ONT methylation
  calls for genome-browser visualization. User will supply their own existing short-read BigWig
  generation scripts next session to adapt rather than building this from scratch.

**Produced:** `notebooks/qc/QC06_pmd_windowed_methylation.ipynb` (replaces DRAFT),
`QC09_haplotype_asymmetry_chain_gaps.ipynb`, `results/qc/data/{qc06_pmd_domains_chr20_pilot,
chain_gaps_raw, chain_gap_summary_all_gaps, chain_gap_summary_sv_ge50bp,
chain_gap_summary_large_ge1kb, postliftover_hap_balance}.tsv`,
`results/qc/figures/qc06_pmd_windowed_methylation_tracks.png`.

**Open / next:**
1. **QC07 still running** — the direct input to per-ancestry CpG-gain/loss and corrected
   CpG-disrupting-SNP rate; check on it first when resuming.
2. **Post-liftover balance check needs AFR/SAS coverage** once P03 progresses further — the
   current "no connection between chain-gap asymmetry and CpG imbalance" finding is untested in
   exactly the population where it matters most.
3. **HMM PMD caller + BigWig generation** — next session, once user supplies their existing
   short-read scripts to adapt.
4. Rare/private-variant extension (per-donor private CpG/SNP fraction, corrected common-vs-rare
   het-site classification via gnomAD) — still planned, not yet built, blocked on QC07/QC08
   scope and the gnomAD query approach above.
5. Carried over: `scripts/github/sync_to_github.sh` still doesn't cover `scripts/download/` or
   `scripts/harmonize_hg38/all_donors/`.

**If resuming, read:** this entry, then check `qstat -u terencew` for QC07/G01-G03 status, then
`notebooks/qc/QC06_pmd_windowed_methylation.ipynb` and `QC09_haplotype_asymmetry_chain_gaps.ipynb`
for the two real findings from this entry.

---

## 2026-09-16 — same day, still later: real HMM PMD caller ported from user's short-read pipeline (revises QC06 substantially), variance-decomposition and bigWig work launched, session-persistence risk identified

**State at start:** the entry above had QC06's threshold-based PMD check and QC09's 15-donor
chain-gap pilot. User asked to scale QC09 to all 202 donors, port their own short-read WGBS
PMD-calling pipeline (`scripts/call_pmds/wgbs_template/`) to this project's modbed data, and
(later in the same session) build a variance-decomposition model and a bigWig export pipeline.

**Decisions made / findings:**
- **Real HMM PMD caller ported and run — substantially revises QC06's earlier finding.**
  New pipeline at `scripts/call_pmds/all_donors/` (`A01a_modbed_to_methcounts.sh` →
  `A01c_pmds.sh`, using `dnmtools sym`/`dnmtools pmd`, jobs **14771403** (methcounts, 404 tasks =
  202 donors × 2 haps) → **14771405** (PMD calls, held on the former)). Validated on the same 3
  donors QC06 used (chr20): **NA19338** (global meth 0.523, low-tail) shows real PMDs (87
  domains, 42.0Mb, 64.2% of chr20 — broadly consistent with QC06's 78.3% threshold estimate) but
  **HG01981** (median, 0.648) and **HG04187** (high, 0.723) show **zero** real PMDs, vs. QC06's
  threshold-based 69.1%/22.4%. **QC06's threshold heuristic was too permissive** — the real HMM
  finds PMDs concentrated in the extreme low-tail donor only, not a smooth gradient across the
  cohort. QC06's conclusion should be treated as superseded pending a full write-up reconciling
  the two.
  - Strand/symmetrization resolved empirically, not assumed: modbed offset sign is fixed by read
    alignment strand (verified — every `-`-strand read has all-negative offsets). Initially
    assumed no merge step was needed from a surface reading of the template (its allc source was
    already pre-merged upstream by methylpy, `*.CGN-Merge.allc.tsv.gz` — that precedent doesn't
    transfer to modbed) — caught by actually reading `dnmtools pmd -help`, which requires
    strand-collapsed input, and empirically confirmed the +/- pairing (3803/4255 plus-strand
    positions had a matching minus-strand call at position+1 in a real 5kb test region).
  - Deliberately used raw (unfiltered) modbed/methcounts, NOT P03's het-filtered hg38 matrix, and
    stayed in native per-haplotype assembly coordinates (matching QC06) rather than lifting —
    the `k≥1` filter would manufacture PMD-like coverage gaps from filtering artifacts, not real
    biology.
  - Found and worked around a `dnmtools` bug: `-S`/`-r`/`-p` (summary/posteriors/params outputs)
    throw an error regardless of argument order on this build; dropped, `-o`'s PMD bed is the
    only output used.
- **Global methylation distribution re-examined for "how many donors have PMDs"**: only 1/219
  donors sits below 0.55 (NA19338 itself), 21/219 below 0.60, bulk (187/219) in 0.60-0.70.
  Working hypothesis, not yet confirmed: PMDs may concentrate in a small low-tail subset (~21
  candidates) rather than being spread continuously — pending the full 202-donor array.
- **Age/passage vs. low-methylation donors, checked directly**: age remains unanswerable (no
  data anywhere). Passage: **11/21 low-methylation donors are missing passage data entirely,
  including NA19338 itself** (the one confirmed real PMD-positive donor). Among the 10 with
  data, mean passage is *lower* (4.60) than the rest of the cohort (4.98), and the full-cohort
  correlation (n=157) is weak and wrong-signed for the "more passage → more PMD-driven
  hypomethylation" story (r=0.196, positive, R²=0.038) — consistent with the earlier finding that
  passage's apparent effect was likely an ancestry-confound artifact, not a real independent
  driver. Whatever drives the low-tail/PMD donors, it isn't simply culture duration.
- **User supplied the correct `igvf_pgp` reference path** (`reference/igvf_pgp/start_merged.filt_mcg.sorted.bed.gz`
  — their own merged fibroblast PMD set) after an earlier session guess at the path was wrong and
  confirmed not to exist. Launched a fork to compare NA19338's PMDs against it (expects modest
  overlap, LCL vs. fibroblast, but wants it measured), plus residual-methylation-in-other-donors
  and standard PMD QC (gene enrichment via `gencode.v43.autosome.dedup.names.bed.gz`, size
  distribution, genomic coverage) — still running as of this entry.
- **Variance-decomposition analysis launched**: linear model of cross-donor methylation variance
  (500bp bins, hg38-lifted) against CpG-island overlap (`cpgIslandExt.hg38.bed`), gene-body
  overlap, and PMD status (donor-specific, represented as fraction-of-donors-in-PMD per bin) —
  testing the CpG-island-low-variance / PMD-high-variance hypothesis directly. Still running.
- **BigWig export pipeline launched**: reasoned through three candidate data sources (raw
  modbed, the new methcounts files, P03's filtered hg38 matrix) and is using unfiltered
  methcounts lifted to hg38 via existing chains — same reasoning as the PMD caller, avoiding the
  het-filter artifact. Still running/validating on the same 3 pilot donors as of this entry.
- **HPRC2 Hi-C/RNA/Fiber-seq download URLs** (user's future validation interest, explicitly
  low-priority): confirmed real per-sample file names via direct bucket listing (not assumed) —
  `hg38.hic` (Hi-C), `expression.{plus,minus}.hg38.bw` (RNA/Kinnex), `fiberseq.PacBio.hap{1,2}.modbed.gz`
  (chromatin accessibility — **not ATAC-seq, HPRC2 doesn't offer ATAC**, only Fiber-seq, and only
  for 21/229 samples). Wrote `tsv/meta/hprc2_hic_rna_fiberseq_urls.tsv` (229 rows, pattern-
  constructed like `A01a`'s existing URLs — not per-sample verified beyond HG00097, some will
  404 for samples missing that assay).
- **Session-persistence risk identified**: this session's interactive shell is a `QRLOGIN` job
  (14764647, `h_rt=24:00:00`, started 05:57:45) — **expires ~05:57 tomorrow**. Several
  long-running fork computations (QC09's full-202-donor scale-up, parts of the variance-decomp
  and bigWig work) are running as plain background bash processes in this shell, NOT via `qsub`
  — they will die if this session expires before they finish. Jobs actually submitted via `qsub`
  (the PMD array 14771403/14771405, dipcall 14766207/10/13) are safe regardless. Worth converting
  any interactive-shell long-runner to `qsub` if it won't finish well within the remaining window.
- **LCL/PMD literature, restated for the record** (already found earlier this session, restated
  here since asked again): Frontiers in Genetics (10.3389/fgene.2017.00076) and corroborating
  sources confirm EBV-transformed LCLs specifically (not just cancer generally) reliably produce
  PMDs — this cohort should be treated as PMD-positive by strong prior, now also confirmed
  directly in NA19338. Plagnol et al. 2008 (PMC2494943) separately established LCL clonality
  (≥22% effectively monoclonal) as a distinct confound. No dedicated large-scale curated
  cross-cell-type PMD atlas was directly verified this session (the user's own fibroblast set at
  `reference/igvf_pgp` is the concrete curated comparison point in hand) — if a broader curated
  resource (e.g. a multi-methylome PMD survey) is wanted, that needs a dedicated literature check,
  not assumed from memory.

**Produced:** `scripts/call_pmds/all_donors/{A01a_modbed_to_methcounts,A01c_pmds}.sh` (+ arrays),
`tsv/meta/hprc2_hic_rna_fiberseq_urls.tsv`, `results/qc/data/het_site_density_by_superpop.tsv`
(from the prior entry, referenced here). Jobs: 14771403/14771405 (PMD calling, running),
QC09 full-scale rerun (in-shell background process, not qsub — at risk, see above), variance-
decomposition and bigWig forks (running).

**Open / next:**
1. **Reconcile QC06 vs. the real HMM PMD calls** — QC06's threshold-based conclusion is now
   superseded for the median/high donors; needs an explicit update or a follow-up note once the
   full 202-donor HMM run lands.
2. **Convert at-risk in-shell background processes to `qsub`** if they won't finish inside the
   ~10.5-hour remaining session window (as of this entry) — real risk of losing partial progress.
3. Fibroblast-overlap, variance-decomposition, and bigWig forks all still running — check and
   fold results in next session if not done by session end.
4. Carried over: rare/private-variant extension via gnomAD (still not built), `sync_to_github.sh`
   still doesn't cover `scripts/download/`/`scripts/harmonize_hg38/all_donors/`.

**If resuming, read:** this entry first — it has the session-persistence risk and the QC06-vs-
real-PMD-caller correction, both load-bearing for trusting anything else in this project right
now.

---

## 2026-09-15 — haplotype-assignment confound investigated end-to-end: real, but not fixable from modbed coordinates alone; het-count filter built, validated, and scaled genome-wide

**State at start:** project had 202/229 assemblies + harmonized modbeds downloaded (from the
2026-09-08 entries) but no analysis pipeline. `md/20260915_haplotype_assignment_modbed_fix.md`
raised a specific concern: HPRC2's hap1/hap2 modbeds are produced by mapping ONT reads to
haplotype-specific assemblies, but for allele-specific methylation (ASM) it only makes biological
sense to trust a read's hap assignment if that read actually spans a site where the two haplotype
assemblies differ (a het site) — a read that maps to hap1 tie-broken by chance carries no real
phasing information. The doc proposed filtering modbed reads by a het-site count threshold `k`.

**Decisions made:**
- **Confirmed the premise is real**, not a hypothetical — orthogonally verified via a forked
  investigation plus the HPRC2 paper's own Methods text: reads are assigned to whichever
  haplotype assembly they align to best, with ties broken arbitrarily, not by genotype-informed
  read-level phasing.
- Built the chr15 pilot as the test case (`scripts/harmonize_hg38/pilot/`), validating against
  98 chr15 imprinted DMRs from Zink et al. 2018 (`data/zink_2018_supp5_pofo_dmrs.csv`) as a
  positive control — these are loci where hap1-vs-hap2 methylation divergence is expected to be
  large and real, so a working filter should raise measured purity/separation on this set without
  needing new ground truth.
- **H01** (`H01_call_hap_het_sites.py`): calls the actual hap1-vs-hap2 het sites per donor/chrom
  via direct assembly-vs-assembly alignment (`minimap2 -cx asm5 --cs` + `paftools.js call`) —
  answers the user's own question directly: het sites come from the assemblies, not the BAMs,
  since modbeds only carry base-mod calls, not genotype. Fixed a bug where `paftools.js call`'s
  mixed `R`(region)/`V`(variant) output lines were treated as uniform BED3 — filtered to `V` only.
- **H02** (`H02_filtered_locus_matrix.py`): filters modbed reads by (a) het-site count ≥ `k`
  spanned by that read and (b) HMMFlagger assembly-reliability overlap (excludes reads over
  Col/Dup/Err/NNN-flagged assembly regions). Fixed a silent-zero bug — `load_het_sites` kept the
  modbed's full `SAMPLE#hap#accession` contig name while the chainmap/lookup used bare accessions,
  so every read's het count silently came back 0 (100% filtered, indistinguishable from "k is too
  strict" without the fix). Also fixed a severe perf bug (ChainMap rebuilt per-call instead of
  cached — 25+ min for 980 calls dropped to 26s once cached).
- **H03** (`H03_validate_k_sweep.py`): swept `k` from 0 to 30 (`[0,1,2,3,5,10,20,30]`) against the
  98 chr15 DMRs, using a purity metric (`max(n_high, n_low) / n_reads` at a 0.5 methylation-
  fraction split, `MIN_READS=3`). **Finding: essentially no benefit from raising k** — max
  locus-level purity change across the entire k range was only 0.08 across 196 (locus, hap)
  pairs. The filter removes reads that can't be trusted, but doesn't resolve the residual
  disagreement.
- **H04** (`H04_epiallele_classification.py`, new): asked *why* purity plateaus — used GMM(k=1)
  vs GMM(k=2) BIC comparison per (region, hap) to classify each locus as genuinely bimodal
  (real ASM signal, purity metric applies) vs. unimodal (not actually imprinted/allele-specific
  at this locus, purity metric is meaningless there). **Only 32.7% of test loci are genuinely
  bimodal.** Even restricting to those, per-read disagreement persists — this is a real biological
  ceiling (equivalent to an "epiallele pattern" question, not a filtering-parameter question), not
  fixable from modbed coordinates because the missing 19-25% of resolving power requires read
  sequence/genotype information the modbed format doesn't carry.
- **Considered and rejected** an alternative: re-deriving everything from raw unaligned HiFi BAMs
  (`tsv/meta/igsr_HPRC2.tsv`) using `asm_lr`'s original from-scratch methodology (align + call
  methylation from kinetics + WhatsHap haplotag), scaled to 230 donors. Rejected as much more
  expensive for a problem that's a real ceiling, not a pipeline bug — the marginal few % of
  reads it might recover isn't worth re-running the full stack at 230-donor scale.
- **Decision on k**: default to `k=1` (any het-site-spanning requirement at all), not the doc's
  suggested k≥2-3 — since H03 showed no measured benefit from going higher, and lower k retains
  more reads/coverage for the same purity.
- **G01** (`G01_call_hap_vs_hg38.py`, new, separate track): implemented hap-vs-hg38 variant
  calling (dipcall-style) for the pilot, since projecting hap-specific coordinates to hg38 is
  independently useful for this project beyond the filtering question. `dipcall-aux.js vcfpair`
  turned out to require GT:AD format that a pure assembly-vs-reference comparison doesn't produce
  (no read depth concept applies) — replaced it with a custom `merge_diploid()` that phases the
  hap1 and hap2 VCFs directly into one diploid VCF, dropping REF-mismatch positions (indel-
  representation disagreements between the two haplotype calls, 1650 dropped in the pilot test).
- **P03** (`P03_build_pilot_matrix.py`, rewritten): now calls H02's filtering logic directly
  instead of orchestrating raw P01+P02 output (H02 is a strict superset of what P01+P02 did).
  Fixed a correctness bug: a (sample, chrom) missing H01 het-site output must be **skipped
  entirely**, not silently treated as "zero het sites everywhere" (which would filter out 100% of
  that locus's reads and look like real "0 kept" data rather than a missing-input gap). Original
  unfiltered 229-donor P01/P02 output preserved untouched at `results/pilot/per_sample/`; new
  filtered output goes to a separate `results/pilot_filtered_k1/per_sample/`.
- **Scaled H01 genome-wide**: submitted `scripts/qsub/H01_call_hap_het_sites_array.sh` as
  **job 14756945**, 4,444 tasks (202 resolved-assembly donors × chr1-22), `-tc 60` throttle,
  `h_data=8G,h_rt=1:00:00` per task. This is the input P03/H02 need to scale past the single
  HG00097/chr15 case they've been validated on so far. Still running as of this entry (51/4444
  tasks concurrently active under the throttle, rest queued).

**Produced:** `scripts/harmonize_hg38/pilot/{chainmap.py, H01..H04, G01, P01..P04}.py`,
`scripts/qsub/{H01_call_hap_het_sites.sh, H01_call_hap_het_sites_array.sh,
G01_call_hap_vs_hg38.sh, P03_build_pilot_matrix.sh}`, `txt/samples/h01_sample_chrom_tasks.txt`
(4,444-row task list), job 14756945 (in progress).

**Open / next:**
1. **Monitor job 14756945** — once it finishes, P03/H02 can scale from the single HG00097/chr15
   case to the full 202-donor × 22-autosome grid. This is the direct prerequisite for building a
   real donors × CpGs ASM matrix at scale.
2. Build out the DMR-level ASM analysis over the full panel once the filtered matrix exists —
   this is the actual paper-facing deliverable; H01-H04 were all pilot/validation work to decide
   *whether* and *how* to filter, not the analysis itself.
3. **`scripts/github/sync_to_github.sh` gap**: only rsyncs `docs/`, `data/`, `scripts/qsub/` —
   does not sync `scripts/harmonize_hg38/` (all the H01-H04/G01/P01-P04 pilot code lives there).
   Should be fixed to include it (excluding the heavy `tmp_chr*/` intermediate dirs, ~24GB each,
   which must stay out of the git mirror) so the mirror doesn't miss the actual pipeline code.
4. Write up the null-ish k-sweep + epiallele-bimodality result plainly in `README.md`/docs — it's
   a real, defensible finding ("k=1 is sufficient; residual disagreement is a genuine biological
   ceiling, not a filtering gap") not a dead end, and should be framed that way when this becomes
   part of the paper's methods/limitations section.

**If resuming, read:** this entry, then `md/20260915_haplotype_assignment_modbed_fix.md` for the
original question, then `scripts/harmonize_hg38/pilot/H02_filtered_locus_matrix.py` and
`H04_epiallele_classification.py` for the two pieces of code that actually answer it.

---

## 2026-09-08 — project founded; HPRC2 discovered, proposal drafted, sample manifest built

**State at start:** `asm_lr` (the predecessor project) had, the day before, locked a paper scope
that explicitly *cut* ancestry-stratified population-genetics claims as underpowered (6 sane
donors, 2 AFR). No plan existed to revisit that.

**Decisions made:**
- HPRC2's public epigenome resource (`hprc-epigenome` S3 bucket, `epigenome.humanpangenome.org`)
  makes the cut population-genetics angle viable again — 229 usable donors vs. 6, AFR now the
  largest superpop group (63) rather than the smallest (2).
- This is additive, not a restart: all 18/18 and 30/30 of `asm_lr`'s donors are present in the
  229. The ASM caller and its calibration fixes, and the validation methodology (imprinted-DMR +
  meQTL enrichment as positive controls), carry over directly.
- Checked HPRC data-use/publication policy: no embargo, no scooping restriction, no pre-review —
  only attribution, an open-access/preprint expectation, and NASEM-guided care around
  superpopulation-label framing.
- Scope stays ASM-only for now — HPRC2 also offers Hi-C, Iso-Seq, Fiber-seq per donor; explicitly
  not pursuing those to avoid becoming an everything-pangenome paper (see README "Explicitly out
  of scope").

**Produced:** `README.md` (proposal), `docs/data_sources.md` (how every S3 path/format claim was
verified), `data/hprc2_sample_manifest.tsv` (229 samples × population/superpop/file
sizes/assembly URLs/overlap flags with `asm_lr`'s cohorts).

**Open / next:**
1. modbed file's actual column/offset encoding doesn't match the published spec (Zhou et al. 2023)
   on a real example row — needs `modbedtools` source, not guessing (`docs/data_sources.md` §5).
2. 30/229 samples' assembly S3 path unresolved (naming-pattern mismatch); 8/229 samples' ancestry
   unmapped against the 1000G panel — both listed by ID in `docs/data_sources.md`.
3. **Not yet checked: whether HPRC2's ONT data is basecaller/chemistry-homogeneous across all
   229.** This is the exact confound that forced `asm_lr` to exclude 12 of its original 30
   donors. If HPRC2 harmonized re-basecalling, that exclusion (and possibly the phase-1 HG01258
   exclusion) may no longer apply — but this must be verified, not assumed.
4. No pipeline code written yet — this session was data discovery and proposal-writing only.

**If resuming, read:** `README.md` in full, then `docs/data_sources.md` §4–5 before writing any
download or parsing script.

---

## 2026-09-08 — same day, later: two open items resolved, download script written, layout fixed

**State at start:** the entry above shipped a proposal with two explicit open questions (modbed
offset encoding, ONT chemistry homogeneity) and 199/229 assemblies resolved. Also: the docs/data
above had been written directly into `github/asm_lr_hprc2/` — the exact working-dir-vs-mirror
divergence pattern flagged as a problem in that same entry, just repeated immediately.

**Decisions made:**
- **modbed format resolved.** Pulled `modbedtools` source (`github.com/lidaof/modbedtools`,
  not just its README) and read the offset-construction code directly, since an earlier
  `WebFetch`-based paraphrase of the format guessed the negative offsets were a "confidence
  score" — wrong. They encode strand/stored-base of each call (`docs/data_sources.md` §5).
  Lesson: for a claim that will drive parser code, read the source, not a summarized paraphrase
  of docs about the source.
- **ONT chemistry homogeneity — largely resolved.** Checked raw-data folders for all 229 samples
  directly (not the paper, which was unreachable — PDF text extraction failed, HTML fetch
  429'd twice). Result: 204/229 have a harmonized Dorado `sup5.0.0` re-basecall, including
  **100% of both the formerly-Guppy (146) and formerly-Dorado-0.6 (56) samples**. Strong evidence
  the chemistry confound that cost `asm_lr` 12 donors is resolved for ~89% of this cohort.
  Real caveat kept in the docs, not dropped: haven't confirmed the epigenome modbeds were
  generated *from* this re-basecall, and a shared model doesn't fully erase an R9.4.1-vs-R10.4.1
  pore-signal difference even if it is the source.
- **Assembly resolution: 199 → 202/229.** The 3 gained were a version-suffix issue (`v1.1.0` not
  `v1.0.1`) — fixed by pattern-matching the directory listing instead of guessing more version
  strings. The remaining 27 have no release2 assembly file at all (not a naming issue), and
  overlap 25/27 with the samples missing the harmonized basecall — treated as one shared
  "not yet fully processed" stratum rather than two separate exclusion lists.
- **Fixed the working-dir/mirror layout.** Moved `README.md`, `JOURNAL.md`, `docs/`, `data/` out
  of `github/asm_lr_hprc2/` into the working directory (`project_ideas/asm_lr_hprc2/`), which is
  now canonical — matching how `asm_lr` itself is laid out (working dir = compute-facing +
  real content, git mirror = derived curated copy). Wrote `scripts/github/sync_to_github.sh`: rsyncs
  the curated set (README, JOURNAL, docs/, data/, scripts/qsub/) into the mirror and stages a
  commit (never auto-commits without a message, never pushes). This is the concrete answer to
  "how do we stop the mirror and working dir from diverging" — one direction of flow, one
  script, run it before ending a session rather than hand-editing the mirror.
- Wrote `scripts/qsub/A01a_download_hprc2_files.sh` (+ `test/` variant, `ID=1` hardcoded,
  matching `asm_lr/scripts/download/D01`'s convention) — one SGE array task per manifest row,
  downloads hap1/hap2 ONT modbed+index and, where resolved, hap1/hap2 assembly FASTA with an
  opportunistic md5 check. Skip-if-exists, writes to a `.partial` path first so a killed task
  never leaves a file that looks complete but isn't.

**Produced:** `docs/data_sources.md` §5 and §7 rewritten with resolved findings;
`data/hprc2_sample_manifest.tsv` regenerated (adds `has_harmonized_sup5_basecall`,
`had_guppy_raw`, `had_dorado06_raw` columns; assembly URLs updated to 202/229 resolved);
`scripts/qsub/A01a_download_hprc2_files.sh` + `test/` variant; `scripts/github/sync_to_github.sh`.
Ran the test variant (task 1, `HG00097`) live to validate the script end-to-end rather than
just eyeballing it.

**Open / next:**
1. Confirm whether `hprc-epigenome` modbeds were generated from the `sup5.0.0` re-basecall
   (provenance not yet directly established, see `docs/data_sources.md` §7).
2. 8/229 samples still unmapped to a superpopulation (not in the 1000G 3202-panel).
3. 27/229 samples (no assembly + mostly no harmonized basecall) — decide whether to exclude as a
   documented stratum (recommended) or investigate further before the full download runs.
4. Once the test download finishes and is spot-checked, submit the full `qsub -t 1-229` array.

**If resuming, read:** this entry, then `docs/data_sources.md` §5 and §7.

---

## 2026-09-08 — same day, session close: full download array submitted, pushed to GitHub

**State at start:** download script validated on 1 sample (task 1, `HG00097`, live test — 2.2GB
modbeds ×2 + ~880MB assemblies ×2, md5s passed). User reviewed the storage estimate (~1.4TB for
all 229) and approved submitting the full array as-is, into the project directory (not scratch).

**Decisions made:**
- Checked quota before submitting: `/u/project/cluo` is at 95.4% (572.5TB/600TB), but ~27TB
  headroom is far more than the ~1.4TB this job needs — safe to proceed.
- Submitted `qsub scripts/qsub/A01a_download_hprc2_files.sh` → **job 14708362, 229 tasks**.
  Task 1 will just skip (files already exist from the earlier test) — idempotent by design, no
  wasted re-download.
- Committed and pushed the working directory's curated set to
  `git@github.com:terencewtli/asm_lr_hprc2.git` (remote already existed, pre-created) via
  `scripts/github/sync_to_github.sh`.

**Produced:** first commit to `github/asm_lr_hprc2`, pushed. Contents: `README.md` (proposal),
`JOURNAL.md` (this file), `docs/data_sources.md`, `data/hprc2_sample_manifest.tsv`,
`scripts/qsub/A01a_download_hprc2_files.sh` (+ `test/` variant). Data itself (`modbed/`,
`assemblies/`, `logs/`) intentionally stays out of the mirror, per the "trimmed mirror"
convention `asm_lr` already uses.

**Open / next (carried over, nothing resolved this entry beyond submitting the job):**
1. **Check `qstat -u terencew` and `logs/A01a_download_hprc2_files.14708362.*` next session** —
   229 tasks downloading ~1.4TB will take a while; some tasks may fail (network, transient S3
   errors) and need re-submission (`qsub -t <failed task IDs> ...`) — the script is idempotent
   so re-running the whole array is also safe/cheap if easier than picking out failures.
2. Confirm whether `hprc-epigenome` modbeds were generated from the harmonized `sup5.0.0`
   re-basecall (provenance still not directly established).
3. 8/229 samples still unmapped to a superpopulation; 27/229 samples (no assembly + mostly no
   harmonized basecall) still need a documented-exclusion decision.
4. No parsing/calling pipeline written yet — this project is still at the data-staging stage.
   Next real step is probably a modbed→per-CpG-call parser (format now understood, §5) once the
   download array finishes.

**If resuming, read:** this entry first (job status!), then `README.md` for full context.
