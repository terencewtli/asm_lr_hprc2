# PROGRESS — informal running log

Operational companion to `RESULTS.md` (findings) and `JOURNAL.md` (chronology). This file is the nitty-gritty:
what has actually been produced, how many of x/y tasks finished, what is missing or stale, and
what a future session must fix before trusting a directory. Update it whenever jobs land; keep
it terse. `RESULTS.md`, `JOURNAL.md`, `JOURNAL.archive.md` and this file live only in this git
mirror — edit them here; `sync_to_github.sh` does not copy them from the working directory.

Last updated: 2026-09-18 ~19:30

## Output inventory

| what | path | count | status |
|---|---|---|---|
| CpG methcounts (native) | `data/pmds/<s>/<s>_hap<h>.cpg.methcounts.tsv.gz` | 404/404 | complete |
| Genome-wide PMD calls (native) | `data/pmds/<s>/<s>_hap<h>.pmd.bed` | 404/404 | complete |
| Per-CpG hg38 bigWigs (meth + depth) | `data/pmds/<s>/<s>_hap<h>.hg38.{meth,depth}.bw` | 402/404 | complete (HG00272 has no chain) |
| 10kb bins + hg38 PMD intervals | `results/meth_bins/per_hap/` | 402 | complete |
| Bin matrix / PCA / domain frequency | `results/meth_bins/` | — | complete (B01b) |
| LCL vs fibroblast boundaries | `results/meth_bins/boundaries/` | — | complete (B02a) |
| PMD QC + mQTL + variance | `results/meth_bins/qc_genetics/` | — | complete (B03a) |
| RT / LAD annotation | `results/meth_bins/annotations/rt_lad_10kb.tsv.gz` | 263,774 bins | complete (B06a) |
| PMD metagene profiles | `results/pmd_metagene/per_hap/` | 402/402 | complete (A01f) |
| Variant density per 10kb | `results/meth_bins/variant_density/` | 22/22 | complete |
| ASM calls per (sample, chrom) | `results/asm/calls/` | 4420/4444 | complete except: 22 HG00272 (no chain, permanent) + **HG02280 chr13/chr16** — `data/het_snps/tmp_chr{13,16}/HG02280/` has the hap FASTAs but no `*_vs_*.paf` / `het_in_hap1.bed`, so H01 died mid-task there; rerun H01 for those two chroms, then C02a tasks 1949/1952 |
| Het-filtered per-CpG counts | `results/asm/cpg/` | 4413 | complete |
| **Donor × CpG matrices** | `results/asm/cpg_matrix/<chrom>.{all,hetfilt}.npz` | 44/44 | complete (C03a) |
| Genome-wide ASM merge | `results/asm/genome/` | 201 donors | complete (C04a) |
| ASM empirical null (chr20) | `results/asm/null/` | 199/202 | complete; summary `results/qc/data/asm_null_vs_real_chr20.tsv` |
| XIST promoter skew | `results/qc/data/xist_promoter_skew.tsv` | 96 females | complete (M03) |
| RNA markers + EBV | `results/qc/data/rna_markers_wide.tsv` | 200 donors | complete (R01a); 29 donors have no RNA file |
| Genome-wide per-donor VCFs | `data/vcf/per_donor_gw/` | 201/202 | complete; cohort merge done: `data/vcf/cohort_gw/all_donors.snps.vcf.gz` (870 MB + .tbi, 2026-09-17 15:25) |
| Solo-WCGW per hap | `results/meth_bins/solo_wcgw/per_hap/` | 402/402 | complete |
| Molecule-level QC (chr20) | `results/qc/data/molecule_qc/` | 201/202 | complete |
| Full 1000G PCA (3,202 samples) | `reference/1000G/pca/g1k_full/pca_result.eigenvec` | 3,202 | complete |
| Figure S1 QC table + PDFs | `csv/figure_s1/`, `pdf/figure_s1/` | 201 donors x 52 cols, 9 panels | complete; panel D re-rendered with genome-wide het counts (14789860, 2026-09-18 00:13) |
| PMD coverage-downsampling test | `results/pmd_downsample/` | — | complete (A01g) |
| Genome-wide het SNVs per donor | `results/qc/data/het_snp_counts_per_donor_gw.tsv` | 201 | complete (bcftools stats on cohort VCF); replaces the 10-donor `QC05_*` pilot file |
| Methylation variance explained (qc17) | `results/qc/data/qc17/` | 229 donors (global), 201 (PMD metrics) | complete (`qc/02b`) |
| Spatial heterogeneity + solo-WCGW clock (qc18) | `results/qc/data/qc18/` | 201 donors, 263k bins, 621 domains | complete (`pmds/05a`, `pmds/06a`) |
| ASM replication: donor tiers + candidates | `results/asm/replication/{donor_tiers.tsv,candidates.tsv.gz}` | 69 discovery donors, 118,796 cpg + 229 Zink candidates | complete (C07a) |
| ASM replication: re-test per donor | `results/asm/replication/calls/` | 201/201 | complete (C07b, 14788189) |
| ASM replication: genotype context | `results/asm/replication/genotype/` | 22/22 | complete (C07c, 14788192) |
| ASM replication: summary + classes | `results/asm/replication/candidates_replication.tsv.gz` | — | complete (C07d, 14788197, 2026-09-18 01:36); `class_summary.tsv`, `nearest_het_by_call.tsv`, `chromhmm_by_class.tsv` |
| ASM method fixes (C08a-d) | `results/asm/model/` | — | see below |
| C08a merged ASM locus sizes | `results/asm/model/asm_locus_size_*.tsv`, `asm_merged_loci.tsv.gz` | 267,208 loci / 69 donors | complete 2026-09-18 |
| C08b adjusted penetrance | `results/asm/model/penetrance_adjusted.tsv.gz`, `donor_propensity.tsv`, `penetrance_null_sensitivity.tsv` | 119,024 regions x 189 informative donors | complete 2026-09-18 |
| C08c lead-variant permutation | `results/asm/model/lead_perm/<chrom>.lead_perm.tsv.gz` | **9/22** (chr11-16, 19-21; all B=1000) | tasks 1-10 of 14808255 died on the bcftools PATH bug (below), resubmitted **14810798**; result so far in RESULTS §12 — 0/5,376 hits fail the permutation and the 1e-4 rule is ~2.1x over-conservative |
| C09a haplotype background | `results/asm/model/hap_background/<chrom>.hapbg.tsv.gz` | 1/22 (chr21, B=300 pilot) | genome-wide **14810779** queued (B=1000); chr21 will be skipped by skip-if-exists — `rm` it to get the full B |
| C08d imprinted-domain reclass | `results/asm/model/candidates_reclassified.tsv.gz`, `imprinted_domains.bed`, `reclass_summary.tsv` | 84 domains, 165 regions moved | complete 2026-09-18 |

### ASM donor tiers (chr20 λ, see JOURNAL)
Authoritative table is `results/asm/replication/donor_tiers.tsv` (202 rows), which is what C07a-d
actually used: **discovery 69 · replication 105 · inflated 22 · excluded 5**. (RESULTS §7's
"70 donors" was an off-by-one against this file; corrected 2026-09-18.) 3 of the 5 "excluded" are
excluded for having **no null run at all**, not for high λ — NA19338 among them, which is the
deepest-domain donor in the cohort, so that exclusion is not random. Also
`results/qc/data/qc16/asm_donor_qc_chr20.tsv` once QC16 runs.

**λ is not just a cutoff — it is a modellable donor property.** log λ_gc regressed on covariates
(199 donors): within-hap read spread R² 0.480, XIST skew 0.189 (95 females), constitutive-domain
depth 0.125, superpopulation 0.049, **sex 0.001 (p = 0.65)**, passage/coverage/read-N50/EBV all
≤ 0.018. Joint: read_sd + depth = **0.659**, + chemistry 0.673 (p = 0.021), + superpopulation
0.678 (all terms p > 0.08). Prefer λ as a continuous covariate over the tier split.

### ASM replication outputs (C07a-d) — read this before quoting any of it

Complete 2026-09-18 01:36. `results/asm/replication/`. Five limits are structural, not bugs; they
constrain what the tables can be asked. **C08a-d (2026-09-18, same day) addressed four of the
five — status flagged inline below. Prefer `results/asm/model/` over the raw C07d columns.**

1. **[ADDRESSED by C08b] The classifier has a hard penetrance floor of 12.4%.** `P_REP = 1e-3` against `P0 = 0.0454`
   with `n_tested_rep = 105` means a region needs **≥ 13/105 donors** to be called anything other
   than `sporadic`/`private`. Everything below is invisible by construction — that is 90.7% of
   candidates. The paper's question is penetrance, so this is the binding constraint, not a
   detail. `P0` is also estimated from the candidate set *including* true signal, which raises
   the floor further; estimate it from matched non-candidate regions instead.
2. **[PARTLY ADDRESSED by C08b — now ±2.4 points, median CI width 0.049] Penetrance resolution is ±8 points.** At n = 105 and p = 0.2, SE = 0.039. Report coarse
   penetrance bands, not point estimates, and never rank individual loci by penetrance.
3. **[ADDRESSED by C08d — imprinted fraction at penetrance >= 0.5 rises 52.9% -> 70.3%] `genotype_indep` is contaminated with imprinting at its high end** — the `zink_10kb` window
   is far too tight for imprinted domains (SNRPN/PWS, ZDBF2, KCNQ1 span 100 kb–2 Mb). 127 of the
   155 regions at penetrance ≥ 0.5 outside that window are in imprinted domains with direction
   consistency ≈ 0.5. Reclassify on imprinted-domain intervals or on `lead_dir_consistency`.
4. **[CHECKED by C08c — the concern was WRONG; see below] `lead_p` is uncorrected for the SNVs scanned.** `C07c` takes `argmin` over a median of 75
   SNVs per region (mean 83, max 1,219) and `C07d` thresholds at 1e-4 with no permutation and no
   LD correction. Naive upper bound ≈ 986 spurious leads across 118,796 candidates vs 18,225
   passing the lead filter alone (~5%). Needs a permutation null before the 6,053
   `genotype_linked` count is publishable.
   **Permutation result (chr21, B=1000): 0 of 319 regions with lead_p < 1e-4 fail at p >= 0.05.**
   The permuted minimum-p has median 0.057, i.e. LD collapses ~75 nominal tests to a handful and
   1e-4 is ~570x beyond the null median. The Bonferroni estimate above assumed independence and
   was wrong. Genome-wide confirmation pending (14808255).
5. **`genotype_linked`'s direction consistency of 1.00 is circular** (`DIR_MIN = 0.8` is one of
   its selection criteria). The non-circular validation is `imprinting` at 0.571 and
   `zink_control` at 0.588 — neither class is selected on that statistic. Quote those.

Two more to state in any methods section:
- **Tier assignment is confounded with yield** (λ vs `n_asm` Spearman 0.974), so discovery donors
  are by construction the lowest-ASM donors. Any ASM caused by the same donor state that raises λ
  is systematically absent from the candidate set. A sensitivity run discovering on a λ-matched
  random subset of replication donors would bound this.
- **Low-read regions are untestable but stay in the denominator.** `MIN_READS = 3`, but complete
  separation gives p = 2/C(n₁+n₂, n₁), so ~14 reads/hap are needed before a perfectly separated
  locus can clear genome-wide BH. Per-hap depth runs 15.6–41.9×, so this dilutes the 0.21% rate
  unevenly across donors. Report a callable-region rate alongside it.

**[RESOLVED by C08a]** `MAX_SPAN = 1000` in `C01a` caps every region at 1 kb, so the
"818 bp / 7 CpG ASM locus" had to be checked against merged loci. It survives: 93.5% of merged
loci are one tile, median merged span 841 bp, p95 1,762 bp. The cap was not binding. Quote
`results/asm/model/asm_locus_size_pooled.tsv`, not the tile width.

**Still open, and now the biggest one:** C08b's answer depends on where the null for `b_i` is
put — 13.9% / 29.6% / 43.6% of candidates called at null quantile 1.00 / 0.75 / 0.50. Use the
conservative default and always show `penetrance_null_sensitivity.tsv` alongside it. Do NOT
iterate the null-set estimate; it has no fixed point (RESULTS §11).

**Genomic control zeroes 12 donors outright** (`donor_propensity.tsv`): 9 of 13 `inflated` and 3
`excluded` call nothing at any candidate after GC, lowest lambda among them 3.84 — not 12.9 as
RESULTS §7 originally said. Effective cohort for anything GC-based is **189 donors**.

### Donor × CpG matrix format (C03a)
`results/asm/cpg_matrix/<chrom>.<source>.npz`, source = `all` (every read, from the bigWigs) or
`hetfilt` (reads spanning ≥1 het site — matches the ASM calls). Arrays:
`pos` (int64, union of hg38 CpGs seen in any haplotype), `cols` (402 `<sample>_hap<N>` strings),
`n_meth` and `n_total` (uint16, n_cpg × n_cols; `n_total == 0` means not covered).
Whole-genome: 22 autosome files per source. **all** 30,870,511 union CpGs (17.6 GB);
**hetfilt** 30,844,989 (17.5 GB). chr20 example: 831,782 union CpGs × 402 haplotypes, median
depth 28 where covered. No chrX/chrY.

## Notebooks: layout, naming, and cohort scope

Grouped by analysis type under `notebooks/<type>/`, named `01a_`, `01b_`, `02a_` … (number =
task group, letter = step within it). Figure notebooks live in
`notebooks/final_figures/figure_<n|sN>/{python,R}`: python exports CSVs to `csv/figure_*`, R reads
those and writes PDFs to `pdf/figure_*` (same split as the lab's YR2_2023 templates).

| notebook | scope | status |
|---|---|---|
| `qc/01a_pclai_sequencing_covariates` | **228 donors / 456 haps** | rebuilt cohort-wide (was 2 donors) |
| `qc/01b_hprc2_supp_seq_qc` | full cohort (HPRC2 S6) | current |
| `qc/01c_actual_vs_reported_ont_coverage` | **402 haps / 197 donors** | rebuilt cohort-wide (was 12 donors, one window) |
| `qc/02a_global_methylation_covariates` | 219 donors | superseded by `qc/02b` (no chemistry term) |
| `qc/02b_methylation_variance_explained_chemistry` | **229 donors (global mCG) / 201 (PMD metrics)** | new 2026-09-17: marginal, chemistry-partial and joint-model R² for 7 metrics x 11 covariates, R9.4.1-only replication (qc17) |
| `qc/03a_xist_skew_qc` | **96 females, promoter window** | rebuilt on M03 output (was gene-body window) |
| `popgen/01a_popgen_af_ld_hetsite_overlap` | 5 donors (het overlap), chr1-3 | partial; het-overlap arm needs gnomAD |
| `popgen/01b_haplotype_asymmetry_chain_gaps` | 15 donors → 202 rerun | current |
| `pmds/01b_pmd_caller_comparison` | **201 donors, genome-wide** | rebuilt; merges the two chr20 pilot notebooks (HMM vs threshold-free) |
| `pmds/02b_variance_inside_outside_pmds` | genome-wide, 201 donors | current |
| `pmds/03a_pmd_size_overlap_and_metagene` | genome-wide, 200 donors | current |
| `pmds/03b_pmd_expansion_metadata` | genome-wide, 201 donors | current |
| `pmds/04a_mechanism_clock_instability` | genome-wide | complete (QC15) |
| `pmds/05a_solo_wcgw_clock` | 402 haps / 201 donors | new 2026-09-17: hap concordance, covariate table, ratio metric (qc18) |
| `pmds/06a_spatial_heterogeneity` | 201 donors, genome-wide | new 2026-09-17: plan E framings 1, 3, 4, 6 (qc18) |
| `asm/01a_asm_calibration_qc` | 201 donors, chr20 | queued (QC16) |
| `final_figures/figure_s1/{python,R}` | **whole cohort QC** | rendered; panel D fix re-running (job 14789860) |

Retired originals live in `notebooks/old/` (kept for provenance, nothing reads them): the 2-donor
read-length/PCLAI pilot, the chr20 windowed-PMD and HMM-comparison pilots, and the chr20/32-donor
variance-by-annotation notebook — the last being a true duplicate of `pmds/02b`, which does the
same decomposition genome-wide with corrected input.

Figure S1 QC table headline numbers (201 donors, cohort-wide): per-CpG depth 30.2x mean
(15.7-41.6); 26.46M hg38 CpGs per haplotype; global mCG 0.662 (0.415-0.748); **2.42M het SNVs per
donor (1.90-3.18M; genome-wide cohort VCF — the earlier "6.77M (6.16-8.77M)" came from the 10-donor
QC05 pilot file and was wrong for Figure S1 panel D)**; median het-site spacing 218 bp; **59 het sites and 68,747 het-CpG linkages per
molecule**; **94.3% of reads carry >=1 het site and 94.5% survive the full filter**; hap1-vs-hap2
global mCG difference 0.003. By superpopulation, AFR leads on heterozygosity (3.11M het SNVs, 97.2%
read retention) and EAS trails (2.21M, 93.1%) — the same axis that drives ASM power.

Cohort-wide results from the rebuilds (2026-09-17):
- PCLAI ancestry PCs for 228 donors; hap1-vs-hap2 PC1 correlation 0.993 (internal consistency).
- Measured vs HPRC2-reported coverage, 197 donors: median ratio 0.92, but **38 donors below 0.8
  and some near 0.4** (HG01361 reported 140x vs 58x measured) — the pilot flagged 2 such donors;
  at cohort scale this is a systematic gap worth understanding (reads not haplotype-assigned or
  not in the modbed).
- XIST promoter skew, 96 females: median 0.393, 36% clonal-like, 14 near-complete (>0.7);
  AFR donors highest (median 0.559).
- PMD caller comparison, 201 donors: the HMM calls ~129k 10kb bins vs ~55k for a threshold-free
  definition (Jaccard 0.28), and the two agree more in deeper donors (r = 0.95 with depth) — the
  HMM is the more permissive of the two.

## Known-stale / do-not-use

- `results/all_donors/per_sample_chrom/` (P03, 166G) — junk calls, unmerged strands, no read
  identity. Superseded by `results/asm/`. Safe to delete.
- `data/vcf/per_donor/` and `data/vcf/cohort/` — **chr21 only** (G02 ran mid-array). Use the
  `_gw` directories once 14780503/14780504 finish.
- `data/pmds/*/*_hap?.methcounts.tsv.gz` (no `.cpg.`) — pre-fix, includes non-CpG junk rows.
- `data/pmds/*/*_hap?.hg38.bw` (6 files) — liftOver-era, 1bp off on flipped blocks.
- `data/pmds/*/*.native.bedGraph`, `*.hg38.bedGraph`, `*.unmapped`, `*.chain` (699G) — temp files
  from the timed-out liftOver bigWig run.
- `results/qc/data/pmd_pilot_chr20/`, QC06/QC10/QC11 — chr20-only, pre-fix inputs.
- Gradient-hub bigWigs not referenced by the current `trackDb.txt` (earlier donor selections).

Deletion commands for all of the above are in `JOURNAL.md` → "Jobs" section.

## Recurring gotchas (cost time at least once)

- **`lead_pos` on disk is CORRUPT in both `C07c` and `C07d` outputs — do not use it.** Both write
  with `float_format='%.5g'`, which truncates genomic coordinates to five significant figures:
  position 606,320 is stored as `6.0632e+05`, so the column cannot identify a variant. Cost a
  full debugging cycle in C09a (87 of 97 chr21 loci silently dropped as "lead not found").
  **The analyses are unaffected** — C07c held the position in memory — only the on-disk column
  is wrong, so RESULTS §8/§13 numbers stand. C09a re-derives the lead from the genotypes instead
  and reproduces the stored `lead_p` exactly (max |dlog10| = 0 over 85 chr21 loci), which is the
  cross-check that proves the re-derivation right. **Fix properly when C07c is next touched:**
  write integer columns with `%d` or drop `float_format` for them. Any other integer column in
  these files (`lead_n_het`, `n_snv_5kb`, ...) is at risk above 99,999 — check before trusting.
- **`bcftools` is not on the PATH inside the job environment.** C08c tasks 1-10 all died with
  `FileNotFoundError: 'bcftools'`. Import `BCFTOOLS` from `C07c_candidate_genotypes` (it holds
  the absolute path `/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools`) rather than relying on
  a bare name. Fixed in C08c/C09a 2026-09-19; tasks resubmitted as 14810798.

- **C08c chr21 holds a smoke-test file with B=200, not B=1000.** It was written by hand while
  developing the script, and the array's skip-if-exists will therefore leave chr21 at the lower
  permutation count. Remove it before or after 14808255 finishes and re-run task 21:
  `rm results/asm/model/lead_perm/chr21.lead_perm.tsv.gz` (then `qsub -t 21 scripts/asm/C08c_lead_permutation_array.sh`).

- `pyBigWig.stats` needs `exact=True`; zoom-level approximation was off >10⁶-fold on long genes.
- `ont_asm_caller` needs `math.comb` (3.8+); the allcools env is Python 3.7 — shim it.
- `Path.unlink(missing_ok=)` is 3.8+; not available in this env.
- UCSC liftOver is far too slow for per-CpG work (~3 h/hap); use `B01a_hap_bins.map_points`.
- `bedToBigBed` needs lexicographic chrom sort; our reference sets are numeric-sorted.
- Array scripts have their own skip-if-exists check — update it when the Python output path
  changes, or every task silently skips (cost one full G02 rerun).
- `pkill -f <pattern>` also matches the wrapper shell running the command.
- R 4.1.0 links against MKL: a qsub job must `source /u/local/Modules/default/init/modules.sh;
  module load R/4.1.0` before nbconvert, or the IRkernel dies with "Kernel died before replying
  to kernel_info" (the missing lib is `libmkl_gf_lp64.so`).
- `domain_frequency_10kb.tsv.gz` and the variant-density tables both have an `n_donors` column;
  merging them silently produces `n_donors_x/_y` and a later KeyError.
- R on compute nodes, second failure mode (after MKL): the system `/lib64/libstdc++` lacks
  `CXXABI_1.3.9`, which the user-library `fastmap` build needs; IRkernel dies the same way. Prepend
  the allcools env `lib/` to `LD_LIBRARY_PATH` for the R step (done in `S1_figures_run.sh`).
- `results/qc/data/QC05_het_snp_counts_per_donor.tsv` is a **10-donor pilot**, not the cohort —
  use `het_snp_counts_per_donor_gw.tsv`.
- In pandas, `df.gt` is the greater-than method: read a column named `gt` as `df['gt']`.
- `zcat f | head -1` under `set -o pipefail` kills the script (SIGPIPE) — read headers in Python.

## Next-session checklist (updated 2026-09-18 ~19:30; the PMD/validation arm has moved to `hprc2_misc`)

**Priorities as of 2026-09-18 night, in order (supersedes anything below that conflicts):**

1. **`C09a` haplotype background — make-or-break.** The lead SNV is a tag, not necessarily the
   cause. Test whether penetrance is predicted by the local HAPLOTYPE carrying the lead allele
   rather than by the lead allele alone. If carriers of haplotype A show ASM and carriers of
   haplotype B do not at the same tag SNP, that is a cis mechanism for the incomplete penetrance
   in RESULTS §13, and it explains the ancestry heterogeneity for free (haplotype frequencies
   differ by ancestry; tag-allele frequencies need not). This is the reason the cohort is
   ancestrally diverse. Everything else is secondary.
2. **`C09b` permutation for the superpopulation-heterogeneity signal.** Turn 7.6%-vs-5% into a
   defensible locus list. chi2 on small per-superpop cells is anti-conservative; permute
   superpopulation labels within lead-genotype strata.
3. **Finish C08c genome-wide** (arrays 14808255 + 14809578) and fold into RESULTS §12, which
   currently cites chr21 only.
4. **Rebuild the ASM narrative around the penetrance continuum** (RESULTS §13 last paragraph).

**Deprioritised:** further caller calibration. lambda is 66% explained (§10), the empirical null
is characterised and its blind spot documented (§7), and C08a-d closed the structural issues.


**Everything submitted on 2026-09-17 has landed.** Nothing of this repo's is queued
(`qstat -u terencew` shows only interactive sessions). What finished overnight:
C07b (201/201) -> C07c (22/22) -> C07d; Figure S1 panel D re-render; QC16 ASM calibration
notebook; U01c hub rebuild + `pmds/03a`; QC15; C05a (after the `has_mqtl` bool fix).

Open, in priority order:
1. **ASM replication stage 4** (why loci fail to replicate) — JOURNAL §2 "ASM replication
   architecture". C07d's output is in and summarized in JOURNAL §0/§1; build the `asm/02a`
   figures notebook off `class_summary.tsv` / `nearest_het_by_call.tsv` first.
2. **Two missing ASM tasks**: rerun H01 for HG02280 chr13 + chr16 (no PAF / het bed), then C02a
   tasks 1949 and 1952.
3. **PMD manuscript, plan E framings 2/5/7** — per-donor boundary calls. NB the Hi-C compartment
   half of framing 5/7 is being built in `hprc2_misc` (`scripts/validate_pmds/hic`, H01a).
4. Cross-repo: de novo meQTLs (hprc2_misc `clusters_k10`) vs PMDs and vs the ASM
   `genotype_linked` class; that run died at the chr9 phenotype task and needs a rerun.

**Submission rule (2026-09-18):** every array gets `#$ -tc` (see `docs/CONVENTIONS.md`).

