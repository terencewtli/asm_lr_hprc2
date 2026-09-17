# PROGRESS — informal running log

Operational companion to `RESULTS.md` (findings) and `JOURNAL.md` (chronology). This file is the nitty-gritty:
what has actually been produced, how many of x/y tasks finished, what is missing or stale, and
what a future session must fix before trusting a directory. Update it whenever jobs land; keep
it terse. `RESULTS.md`, `JOURNAL.md`, `JOURNAL.archive.md` and this file live only in this git
mirror — edit them here; `sync_to_github.sh` does not copy them from the working directory.

Last updated: 2026-09-17 ~16:00

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
| Solo-WCGW per hap | `results/meth_bins/solo_wcgw/per_hap/` | 217/402 | **running** (B05a, job 14779380) |
| Variant density per 10kb | `results/meth_bins/variant_density/` | 22/22 | complete |
| ASM calls per (sample, chrom) | `results/asm/calls/` | 4410/4444 | 12 real gaps resubmitted (14780506); 22 are HG00272 (no chain, permanent) |
| Het-filtered per-CpG counts | `results/asm/cpg/` | 4413 | complete |
| **Donor × CpG matrices** | `results/asm/cpg_matrix/<chrom>.{all,hetfilt}.npz` | 44/44 | complete (C03a) |
| Genome-wide ASM merge | `results/asm/genome/` | 201 donors | complete (C04a) |
| ASM empirical null (chr20) | `results/asm/null/` | 199/202 | complete; summary `results/qc/data/asm_null_vs_real_chr20.tsv` |
| XIST promoter skew | `results/qc/data/xist_promoter_skew.tsv` | 96 females | complete (M03) |
| RNA markers + EBV | `results/qc/data/rna_markers_wide.tsv` | 200 donors | complete (R01a); 29 donors have no RNA file |
| Genome-wide per-donor VCFs | `data/vcf/per_donor_gw/` | 201/202 | complete; cohort merge (G03) running |
| Molecule-level QC (chr20) | `results/qc/data/molecule_qc/` | 1/202 | **queued** (Q01a, job 14781764) |
| Full 1000G PCA (3,202 samples) | `reference/1000G/pca/g1k_full/` | — | **queued** (job 14781790) |
| PMD coverage-downsampling test | `results/pmd_downsample/` | — | **queued** (A01g, job 14781691) |

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
| `qc/01a_read_length_coverage_pclai` | **2 donors** | pilot placeholder; superseded by figure_s1 |
| `qc/01b_hprc2_supp_seq_qc` | full cohort (HPRC2 S6) | current |
| `qc/01c_actual_vs_reported_ont_coverage` | 12 donors, single window | superseded by figure_s1 (measured depth for all 402 haps) |
| `qc/02a_global_methylation_covariates` | 219 donors | current, but predates chemistry — see RESULTS §1 |
| `qc/03a_xist_skew_qc` | 98 donors, wrong window | superseded by `M03_xist_promoter_skew.py` |
| `popgen/01a_popgen_af_ld_hetsite_overlap` | 5 donors (het overlap), chr1-3 | partial; het-overlap arm needs gnomAD |
| `popgen/01b_haplotype_asymmetry_chain_gaps` | 15 donors → 202 rerun | current |
| `pmds/01a_pmd_windowed_methylation` | chr20, 3 donors | superseded (threshold method) |
| `pmds/01b_pmd_hmm_comparisons` | chr20 pilot | superseded |
| `pmds/02a_variance_decomposition_by_annotation` | chr20, 32 donors, pre-fix input | superseded by `02b` |
| `pmds/02b_variance_inside_outside_pmds` | genome-wide, 201 donors | current |
| `pmds/03a_pmd_size_overlap_and_metagene` | genome-wide, 200 donors | current |
| `pmds/03b_pmd_expansion_metadata` | genome-wide, 201 donors | current |
| `pmds/04a_mechanism_clock_instability` | genome-wide | queued (job 14781796) |
| `asm/01a_asm_calibration_qc` | 201 donors, chr20 | queued (job 14781797) |
| `final_figures/figure_s1/{python,R}` | **whole cohort QC** | queued (job 14781798) |

The four "superseded" notebooks are kept for provenance; nothing downstream reads them. Their
questions are answered cohort-wide by figure_s1, `pmds/02b`–`04a` and `M03`.

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

## Next-session checklist

1. `qstat -u terencew`; then check the counts table above.
2. Held/queued: QC15 (14779515), QC16 (14780472), C05a (14779512), U01c hub rebuild (14778024).
3. Solo-WCGW (14779380) and the ASM empirical null (14780466) are the two inputs blocking the
   PMD-clock and ASM-calibration write-ups.
4. Re-run `U01a`/`U01b` after any new bigWigs, and re-execute QC12 for the metagene section.
