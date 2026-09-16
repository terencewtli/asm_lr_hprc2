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
