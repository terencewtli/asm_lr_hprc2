# Data source identification — 2026-09-08

How the S3 paths, sample list, and format claims in `README.md` and
`data/hprc2_sample_manifest.tsv` were established. `epigenome.humanpangenome.org` itself is a
JS single-page app — its content isn't fetchable by URL, so everything below was derived from
the underlying public S3 buckets directly.

## 1. Locating the bucket

`WebSearch` → [Registry of Open Data on AWS: hprc-epigenome](https://registry.opendata.aws/hprc-epigenome/)
→ public, unsigned-request bucket:

```
aws s3 ls --no-sign-request s3://hprc-epigenome/
```

No `aws` CLI on Hoffman2 (checked: not on `$PATH`). Used the S3 REST `ListObjectsV2` API over
plain HTTPS instead, which works identically for a public bucket without credentials:

```
curl -s "https://hprc-epigenome.s3.amazonaws.com/?list-type=2&delimiter=/&prefix=samples/"
```

paginated via `<NextContinuationToken>` for the full 232-sample, 11,355-object listing
(`/tmp/s3_all_keys.txt` in the session that ran this, not checked into git — regenerate with the
loop below if needed).

## 2. Sample directory structure

`s3://hprc-epigenome/samples/<SAMPLE_ID>/<file>`. Confirmed via the data provider's own tutorial
notebook, [`twlab/open-data-examples/get-to-know-hprc-epigenome.ipynb`](https://github.com/twlab/open-data-examples/blob/main/get-to-know-hprc-epigenome.ipynb)
(fetched via `raw.githubusercontent.com`, not the `github.com` blob view — the latter doesn't
serve notebook cell content through `WebFetch`).

Per-sample files actually observed (using `HG00097` and `HG01258` as examples):
`methylation.ONT.hap{1,2}.modbed.gz(.tbi)`, `methylation.PacBio.hap{1,2}.methylc.gz(.tbi)`,
`fiberseq.PacBio.hap{1,2}.modbed.gz(.tbi)` (21/229 samples only), `hap{1,2}.refbed.gz`,
`hap{1,2}_vs_{chm13,hg38}.gz` (per-donor liftover/alignment to reference — useful for cross-donor
coordinate harmonization, not yet used), `{chm13,hg38}.hic` + boundary/insulation tracks,
`expression.{plus,minus}.*.bw`, `CGI.bed.gz`, `HMMFlagger.{ONT,PacBio}.bed.gz`,
`RepeatMasker.hap{1,2}.bb`.

229/232 samples have both `methylation.ONT.hap1.modbed.gz` and `.hap2.modbed.gz` — the 3 without
(`HG002`, `HG02622`, `NA19087`) were excluded from the manifest.

## 3. Ancestry / superpopulation

Not in the epigenome bucket. Most sample IDs are 1000 Genomes IDs, so joined against the public
1000G 30x panel:

```
curl -s "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"
```

(3,202 samples, columns include `Population`/`Superpopulation`). **8 of the 229 usable samples
aren't in this panel**: `HG005, HG01123, HG02109, HG02486, HG02559, HG03471, HG06807, NA21309`.
Not yet resolved — some (`HG002`/`HG005`) are known GIAB reference samples outside 1000G; the
rest are unexplained and need HPRC's own sample metadata (not yet located) rather than an
assumption.

## 4. Assembly (FASTA) download paths

Not in the `hprc-epigenome` bucket at all — assemblies live in the pre-existing
`human-pangenomics` bucket already used by `asm_lr` (see its `csv/meta/cohort_downloads.tsv`,
which had 30 rows with working `assembly_hap{1,2}_s3` paths for the original cohort — used as the
template pattern).

Verified pattern for release-2, non-trio samples via `curl -I` (HTTP HEAD, no download):

```
https://human-pangenomics.s3.amazonaws.com/working/HPRC/<ID>/assemblies/release2/<ID>_hap1_hprc_r2_v1.0.1.fa.gz
https://human-pangenomics.s3.amazonaws.com/working/HPRC/<ID>/assemblies/release2/<ID>_hap2_hprc_r2_v1.0.1.fa.gz
```

**Trio-phased samples use a different naming scheme** (`_pat_`/`_mat_` instead of `_hap1_`/
`_hap2_`, since parent-of-origin is known unambiguously rather than arbitrarily assigned) — e.g.
`HG04228_mat_hprc_r2_v1.0.1.fa.gz` / `_pat_`. The manifest's `assembly_naming` column records
which pattern resolved (`hap` or `trio`) per sample; this distinction matters for the eventual
pipeline (mat/pat samples have a real parent-of-origin label for free, useful for the imprinting
class specifically) and should be tracked, not discarded.

Checked both patterns (`hap1`/`hap2`, `pat`/`mat`) × version suffix `v1.0.1` for all 229 usable
samples via concurrent `HEAD` requests: 199 resolved immediately. For the 30 that didn't, listed
each one's full `assemblies/` and `assemblies/release2/` prefixes directly (same technique as
§1) rather than guessing more filename variants:

- **3 resolved on a second pass** (`HG01978`, `HG02257`, `HG03516`) — same `mat`/`pat` trio
  naming, but version suffix `v1.1.0` instead of `v1.0.1`. **Lesson: don't hardcode a version
  string in the download script — list the directory and pattern-match `*.fa.gz`, the way the
  fallback resolver here does, since the suffix isn't uniform across the cohort.**
- **27 have no assembly file at all** under `working/HPRC/<ID>/assemblies/release2/` — the
  directory is empty or the `release2/` prefix doesn't exist under `assemblies/` for that sample.
  Not a naming problem: `HG005, HG00733, HG01109, HG01243, HG02055, HG02080, HG02109, HG02145,
  HG02723, HG02818, HG03098, HG03486, HG06807, NA18906, NA18940, NA18943, NA18944, NA18945,
  NA18948, NA18959, NA18960, NA18967, NA18970, NA18982, NA19240, NA20129, NA21309`. Final
  resolution: **202/229**.

This 27-sample list overlaps heavily (25 of 27) with the samples missing the harmonized
`sup5.0.0` ONT re-basecall in §7 below — most likely one underlying stratum ("not yet through
this release's full harmonized processing"), not two independent gaps. Recommend tracking and
excluding it as a single documented set, the same way `asm_lr` documented its 12-donor
chemistry exclusion, rather than as two separate ad hoc filters.

Note `NA18959` is one of the 27 — it's also flagged in `asm_lr/md/20260904_cohort_provenance.md`
§9 as having "no confirmed raw-data directory under HPRC `working/`" for the *old* cohort's
raw-BAM download. Given it now also lacks a release2 assembly and a harmonized re-basecall, this
looks like the same underlying gap (this donor's data simply isn't fully deposited yet), not
three unrelated coincidences.

## 5. modbed file format — RESOLVED

Column layout confirmed via the format's original publication: Zhou et al., *Cell Genomics*
(2023), ["Modbed track: Visualization of modified bases in single-molecule sequencing"](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00299-9)
(Wang lab, WashU — same group hosting the HPRC2 epigenome browser). Per that paper: cols 1–3
chrom/start/end, col 4 read ID, col 5 score, col 6 strand, col 7 comma-separated
methylated-base offsets, col 8 comma-separated unmethylated-base offsets, both "relative to the
start position (col 2)".

A pre-existing local test file
(`/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/tmp/tmp.modbed.gz` — not gzip-
compressed despite the name, plain text) has rows like:

```
HG00097#2#CM094075.1  427  7904  146c7282-...  0  -  -6665,-6677,-6689,...  -6713,-6725,-6737,-6743
```

An initial `WebFetch`-based read of the paper guessed the negative offsets were a "confidence
score" — **that guess was wrong**, caught by going to the actual generating code rather than
trusting a natural-language paraphrase of the paper. Pulled `modbedtools` source directly
(`github.com/lidaof/modbedtools`, `modbed/modbed.py`) and read the offset-construction logic
(lines ~100–131): the sign of each offset is **not noise or a score, it encodes which strand /
stored base the call came from**. Concretely, when a modification call's underlying base is
stored as `T` (the reverse-complement representation of a reference `C`, as happens for calls on
the minus-strand copy of a CpG), the offset is emitted with a negative sign; when it's stored as
the forward base, the offset is positive. This lets a single pair of columns report CpG calls
from either strand at consistent genomic coordinates without a dedicated per-offset strand
column — position is `abs(offset)` from col 2, strand-of-call is `sign(offset)`. A read whose
listed CpGs are consistently all-negative (as in the example above) is simply one whose calls
were all stored in `T`/reverse-complement form, which is unremarkable.

This is enough to write a correct parser: `genomic_position = col2 + abs(offset)` for each
offset in col 7 (methylated) or col 8 (unmethylated); `sign(offset)` need only be tracked if the
per-strand distinction matters downstream (for haplotype-level CpG aggregation across a
diploid-assembly-coordinate contig, it likely doesn't, since strand is already fixed by contig
orientation — worth confirming with one real worked example before assuming, not asserting).

## 7. ONT basecaller/chemistry homogeneity — checked, largely resolved

`asm_lr` excluded 12 of its original 30 donors because their raw ONT data was R9.4.1/Guppy-called
rather than R10.4.1/Dorado — a real confound for pooled methylation calling. Rather than trust
the bioRxiv preprint's methods text (couldn't reliably fetch it — see below), checked the actual
raw-data directory structure directly, since `asm_lr`'s own `cohort_downloads.tsv` already showed
the path shape: `s3://human-pangenomics/working/HPRC/<ID>/raw_data/nanopore/<basecaller_folder>/`.

Listed that prefix (via the same `ListObjectsV2` REST technique as §1, `delimiter=/`) for all 229
usable samples. Result: **204/229 (89%) have a `sup5.0.0` (Dorado super-accuracy v5.0.0) folder
alongside their original basecall folder** — and critically, **100% of the 146 samples whose
original data was `guppy_6` and 100% of the 56 originally `dorado0.6.0`-called samples have it**.
This is strong, systematic (not spot-check) evidence that HPRC2 re-basecalled essentially the
whole cohort with one harmonized model regardless of original chemistry batch.

**What this does and doesn't establish:**
- Does establish: a shared-model re-basecall exists, at scale, for nearly everyone — the
  chemistry split that drove `asm_lr`'s 12-donor exclusion is very likely no longer a hard
  blocker for most of the cohort.
- Does NOT establish: that `hprc-epigenome`'s `methylation.ONT.hap*.modbed.gz` files were
  actually generated from this `sup5.0.0` re-basecall specifically, rather than from the older
  basecall. No provenance metadata for this was found in the S3 listing; would need to inspect
  read IDs or a BAM/modbed header, or find a methods description, to close this gap.
- Does NOT fully retire the R9.4.1-vs-R10.4.1 concern even if the provenance link holds — a
  shared basecaller *model* reduces but doesn't necessarily eliminate a real pore/signal-level
  chemistry difference. Recommend an empirical check (does pooled data show a chemistry-batch
  effect) rather than either assuming it's fine or re-imposing the old exclusion by default.

The 25 samples missing `sup5.0.0` overlap 25/27 with the samples missing a release2 assembly
(§4) — treat as one likely-shared "not yet fully processed in this release" stratum:
`HG005, HG02055, HG02080, HG02109, HG02145, HG02723, HG02818, HG03098, HG03486, HG06807,
NA18906, NA18940, NA18943, NA18944, NA18945, NA18948, NA18959, NA18960, NA18982, NA19240,
NA20129, NA21309` plus a few others — full per-sample folder listing in
`data/hprc2_sample_manifest.tsv` (`has_harmonized_sup5_basecall` column).

Tried to corroborate via the bioRxiv preprint's own methods text directly (would be more
authoritative than inferring from directory structure) — `WebFetch` on the PDF returned only
raw PDF object structure, not text; the PDF was saved locally but `Read` couldn't render it
(`pdftoppm`/`poppler-utils` not installed on Hoffman2); `WebFetch` on the `.full` HTML view
was rate-limited (HTTP 429) both times it was tried. The S3-directory-structure approach above
was used instead since it's the primary data, not secondary reporting of it — but if the paper's
methods section is later accessible, it should still be checked, since it may state this
explicitly rather than requiring inference.

## 8. Policy pages

`humanpangenome.org/data-use/` and `humanpangenome.org/publication-policy/`, fetched directly —
summarized in `README.md`. No further verification needed; these are the consortium's own
authoritative pages.

## What wasn't checked

- Whether `hprc-epigenome`'s modbed files were actually generated from the harmonized
  `sup5.0.0` re-basecall found in §7, vs. inferring it from the existence of that re-basecall in
  the separate raw-data bucket — see §7's caveat.
- Actual byte-level parsing of any modbed file (only file sizes were used, as a coverage proxy).
- ~~The `hprc-epigenome` bucket's `files.html` object~~ — checked; it's just a generic
  client-side S3 browser widget (JS that calls the same `ListObjectsV2` REST API used above), no
  format documentation in it.
