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
