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

### ASM donor tiers (chr20 λ, see JOURNAL)
λ < 1.2: 70 donors (use as-is) · 1.2–2: 69 · 2–3: 37 (genomic control) · > 3: 25 (exclude from
genome-wide discovery) · > 8: NA20762 + 2 (exclude outright). Table:
`results/qc/data/qc16/asm_donor_qc_chr20.tsv` once QC16 runs.

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

