# PROGRESS — informal running log

Operational companion to `JOURNAL.md` (which holds conclusions). This file is the nitty-gritty:
what has actually been produced, how many of x/y tasks finished, what is missing or stale, and
what a future session must fix before trusting a directory. Update it whenever jobs land; keep
it terse. Both files live only in this git mirror.

Last updated: 2026-09-17 ~14:30

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
| Solo-WCGW per hap | `results/meth_bins/solo_wcgw/per_hap/` | 181/402 | **running** (B05a, job 14779380) |
| Variant density per 10kb | `results/meth_bins/variant_density/` | 21/22 | chr22 was node-killed; resubmitted (14780505) |
| ASM calls per (sample, chrom) | `results/asm/calls/` | 4410/4444 | 12 real gaps resubmitted (14780506); 22 are HG00272 (no chain, permanent) |
| Het-filtered per-CpG counts | `results/asm/cpg/` | 4413 | complete |
| **Donor × CpG matrices** | `results/asm/cpg_matrix/<chrom>.{all,hetfilt}.npz` | 44/44 | complete (C03a) |
| Genome-wide ASM merge | `results/asm/genome/` | 201 donors | complete (C04a) |
| ASM empirical null (chr20) | `results/asm/null/` | 4/202 | **running** (C06a, job 14780466) |
| XIST promoter skew | `results/qc/data/xist_promoter_skew.tsv` | 96 females | complete (M03) |
| RNA markers + EBV | `results/qc/data/rna_markers_wide.tsv` | 200 donors | complete (R01a); 29 donors have no RNA file |
| Genome-wide per-donor VCFs | `data/vcf/per_donor_gw/` | 0/202 | **rerunning** (14780503 → G03 14780504); first attempt skipped everything (array skip-check pointed at the old path) |

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
