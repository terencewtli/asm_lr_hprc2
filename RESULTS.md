# RESULTS — what we know, and what shows it

Current validated findings only. No history: how each result was arrived at (including analyses
that were retracted along the way) is in `JOURNAL.md`; what has finished running is in
`PROGRESS.md`.

Cohort: HPRC2 lymphoblastoid cell lines, haplotype-resolved ONT methylation.
**202 donors with assemblies, 402 haplotypes, ~30x per CpG per haplotype**, hg38.
Last updated 2026-09-17.

---

**Two manuscripts.** The PMD/domain findings (§2-5) and the ASM findings (§6-7) are being written
in parallel as separate papers; §1 (resource and its technical limits) is shared by both.

## The questions

1. Is this dataset usable as a population-scale haplotype methylome resource, and what are its
   technical limits?
2. What does variation between LCL methylomes actually look like?
3. What causes it?
4. How do these domains compare to PMDs in another cell type?
5. Is domain state genetic?
6. How much allele-specific methylation is there, and are the calls real?
7. What has to be controlled before ASM can be compared across donors?

---

## 1. The resource is usable; the main technical limit is ONT chemistry

**Coverage.** ~30x per CpG per haplotype (mean/median 30.1/30 for NA19338 hap1 over 32.2M
assembly CpGs, 99.93% of them covered); ~35x physical read depth, mean read span 57 kb.

**Coverage is sufficient per haplotype — verified across all 402, not extrapolated.** Mean
per-CpG depth of the methcounts actually fed to `dnmtools pmd`: median **29.7x**, quartiles
27.5–33.3, **min 15.6x, max 41.9x; none below 15x, only 6 haplotypes below 20x**. `dnmtools`
documents ~10x as its recommendation, so every haplotype is ~1.5–4x above it. Coverage barely
predicts the calls (corr with PMD burden +0.09, with domain contrast +0.10).

**Why we call PMDs per haplotype rather than pooling the two into a "diploid" methylome:**
- coverage is not the binding constraint (above), so pooling buys little;
- the haplotypes agree about where domains are (PMD-bin Jaccard 0.91 per donor; mean |hap1−hap2|
  0.026 inside domains vs 0.019 outside), so pooling would not move the map;
- pooling is not free: the two haplotypes are assembled separately, so it requires projecting
  both into a common reference first, and it discards exactly the haplotype resolution the ASM
  arm of the project needs.
A downsampling test (thin each haplotype to 20/15/10/5x and re-call) is running to state the
coverage-sufficiency point empirically rather than by reference to the tool's recommendation.

**Chemistry is the dominant technical covariate and is NOT harmonized.** The modbeds use each
donor's original basecalls: 156 donors R9.4.1/Guppy, 73 R10.4.1/Dorado (HPRC2 Supp S6). Verified
by read ID: NA19338 and HG00099 modbed reads are 6,000/6,000 from their 2022 R9 runs, although
newer R10 data exists for them.

- R10 donors read ~4 points lower in call-weighted global mCG (R² 0.21–0.26); the effect is
  concentrated inside low-methylation domains (in-domain mCG 0.56 vs 0.61).
- It is PC2 of the 10 kb methylation matrix (R² 0.30).
- It does **not** move domain locations: R9-vs-R10 domain-frequency correlation r = 0.994.
- HPRC2 themselves use ONT chemistry as an mQTL covariate but report no effect size.

**Practical rule:** put chemistry in every cross-donor model and replicate headline results
within R9.4.1 alone (n ≈ 150).

Also verified: coverage matches HPRC2's reported values (`figures/qc/coverage_actual_vs_reported.png`),
and the chr1 centromere gap in browser views is an hg38 limitation (modeled satellite / N-gap),
not missing data — the assemblies carry that sequence at normal coverage.

---

## 2. LCL methylomes differ along one continuous axis: how deep their domains are

**PC1 of the 10 kb matrix explains 82% of between-donor variance**, and tracks domain
hypomethylation (PMD burden R² 0.59), not ancestry (0.11).

**The gradient is smooth, not on/off.** Every donor has domains (genome-wide PMD burden
1.10–1.89 Gb, median 1.45; no donor near zero); burden is unimodal; depth is continuous with a
long deep tail.

![domain depth continuum](figures/pmds/qc12_domain_depth_continuum.png)

*Left: each donor's own PMD contrast vs its global mCG — a clean continuum from NA20762 and
NA19338 (deep) to HG04187 and HG02129 (shallow). Right: mCG in constitutive domain bins vs
never-domain bins; all donors sit below the diagonal, i.e. domains are always the lower
compartment, by a donor-specific amount.*

**Domain locations are shared; only depth varies.** 41% of 10 kb bins are domains in ≥90% of
donors; donor-vs-donor PMD overlap is Jaccard 0.87. The deepest donor (NA19338) is an extreme
version of the common methylome, not a different one: 65% of its PMD bins are the constitutive
set. Even the most methylated donors show the same domains, just shallowly.

![PMD metagene by global mCG](figures/pmds/qc12_pmd_metagene_by_global_mcg.png)

*Metagene over constitutive domains and fibroblast PMDs (40 scaled body bins, 200 kb flanks).
Every donor has the same domain shape and boundaries; the curves differ almost only in depth.*

**Most between-donor variance lives in domains, and the rest of the genome follows.**
Constitutive-domain bins are 41% of bins but carry 60% of between-donor variance (all domain
classes: 90%); never-domain bins are 31% of bins and 9.7% of variance. Outside-domain mCG rises
linearly with in-domain mCG (slope 0.27, R² 0.84, no curvature), i.e. it moves at about a quarter
of the in-domain rate, with no threshold where it becomes "normal".

![inside vs outside domains](figures/pmds/qc14_inside_vs_outside_domains.png)

Outside domains, the state-linked variation is concentrated in intermediately methylated
intergenic sequence (R²_state 0.68) and is absent from gene bodies (0.22).

![where outside-domain variance lives](figures/pmds/qc14_state_variance_outside_domains.png)

---

## 3. The domain map is replication timing; the donor-to-donor depth is still unexplained

**Structural cause — settled.** Using GM12878 (an LCL) Repli-seq and lamin-B1 LADs lifted to
hg38 10 kb bins: domain frequency vs replication timing **Spearman −0.82**. In the latest
replicating decile 99.2% of bins are constitutive domains and 62% are LAD; in the earliest, 0.3%
and 5%. LAD overlap by class: 13% (never) → 61% (constitutive). This is the late-replicating,
lamina-associated compartment, as the literature predicts.

**Depth per donor — not explained by anything measured.** Tested and rejected or weak
(partial R² after chemistry, on domain depth):

| candidate | partial R² | verdict |
|---|---|---|
| passage number | ~0.00 | no (and 144/155 recorded lines are p5) |
| donor age | — | not available (Coriell: "Age: No Data") |
| sex | ~0.00 | no |
| superpopulation | 0.04 | weak; drops out with chemistry in the model |
| line established at Coriell vs external | 0.00 | no |
| banking era (NA vs HG ID) | 0.02 | weak, confounded with chemistry |
| **clonality (XIST promoter skew)** | **0.004 (p = 0.54)** | **no** |
| EBV transcription (chrEBV in Kinnex RNA) | 0.014 (p = 0.10) | no |
| proliferation / plasmablast / naive-B / activation / senescence RNA | ≤0.015 | no |
| coverage, read N50 | ≤0.02 | no |

![PMD metrics vs metadata](figures/pmds/qc13_pmd_depth_by_metadata.png)
![PMD metrics by passage](figures/pmds/qc13_pmd_metrics_by_passage.png)

**Clonality is real but is not the driver.** The XIST promoter assay (fixed from an earlier
gene-body window) shows **35/96 female donors (36%) with skew > 0.5**, several essentially fully
skewed — independently reproducing Plagnol 2008's ≥22% monoclonality estimate. It does not
predict domain depth.

**"Expansion" is one axis, not two.** Depth, fixed-threshold breadth, relative breadth and spread
into non-constitutive bins all correlate at r = 0.97–0.99. Domains deepen in place; boundaries do
not move. The open question is therefore what sets a line's cumulative division history, which
nothing in this dataset measures directly. The solo-WCGW clock (3.39M sites annotated; deepens
NA19338's domains from 36.6% to 15.0% mCG) is the sharpest available handle.

---

## 4. LCL and fibroblast PMDs agree on large domains, not small ones

Against the lab's fibroblast (IGVF PGP) PMD set:

| fibroblast PMD size | n | fraction inside LCL PMDs |
|---|---|---|
| <100 kb | 308 | 0.17 |
| 100–300 kb | 391 | 0.26 |
| 300 kb–1 Mb | 437 | 0.58 |
| 1–3 Mb | 220 | 0.90 |
| >3 Mb | 54 | 0.96 |

![fibroblast/LCL overlap by size](figures/pmds/qc12_fib_lcl_pmd_overlap_by_size.png)

Fibroblast PMDs >1 Mb are domains in ≥95% of donors. The overall Jaccard stays ~0.51 (0.52 with
both sets restricted to ≥300 kb) because **LCLs carry ~400–500 Mb of large domains fibroblasts
lack** — not because small domains disagree. Boundaries shared between the two cell types are
gene-enriched (1.13x, permutation p = 0.01); LCL-specific recurrent boundaries are not (p = 0.17).

---

## 5. Genetic regulation is present in domains; the extra variants there are replication timing

**mQTLs are not excluded from domains.** Among TSS-containing 10 kb bins, adjusted for CpG and
TSS count, HPRC2 promoter mQTLs are mildly *enriched* in domain bins (OR 1.18 constitutive, 1.36
variable, 1.47 rare; p < 1e-4) relative to never-domain bins.

**Sequence divergence in domains is a replication-timing effect.** From 790M assembly-vs-hg38
SNV calls across ~200 donors: constitutive-domain bins carry 13.6 SNVs per donor per 10 kb vs
11.1 in never-domain bins, but the domain coefficient attenuates **82% once replication timing
and LAD are in the model** (indels 99%; SVs show no domain association at all). SNV density rises
monotonically from 11.1 to 14.0 across RT deciles. So there is no PMD-specific instability beyond
late replication.

---

## 6. ASM: ~0.2% of regions per donor, and the recurrent calls are the imprinted genome

Read-level test (per-molecule methylation fractions, Welch/Mann-Whitney) over 2,019,218 hg38 CpG
clusters (median 852 bp, 10 CpGs) covering 26.1M of 27.7M autosomal CpGs.

**Rate, in the 70 well-calibrated donors (see §7):**
- median **4,077 ASM regions per donor = 0.21% of tested regions = 0.23% of tested CpGs**,
  median |Δ| 0.26; ASM regions are 818 bp / 7 CpGs.
- union across those donors: **119,926 distinct regions (5.9% of tested)**.

**Comparison:** deCODE (Nat Genet 2024, 7,179 Icelanders, ONT) report 1.2% of CpG units as
candidates and 0.51% validated as ASM-QTLs. Their number is a cohort union validated against
genotype; ours is per donor with no genotype step, so 0.23% per donor is the same order and the
cohort union is the quantity that would need genotype validation to be comparable.

**Recurrence is the structure that matters:** 66% of union regions are called in exactly one
donor, 26% in 2–5, 1% in >20; **407 regions in >50% of donors, 239 in ≥90%**.

**The ≥90% set is imprinting.** 183/239 (77%) lie within 10 kb of a Zink 2018 imprinted DMR;
they collapse to ~41 clusters on chr15 (108 regions, SNRPN/PWS), chr20 (35, GNAS), chr2 (27),
chr11 (14, H19/IGF2), chr7 (11, MEST/GRB10), chr19 (9). Per donor, 43–45% of the 229 Zink DMRs
are recovered with mean |Δ| 0.25. The caller finds the imprinted genome without being told about
it.

---

## 7. Per-donor calibration is mandatory before any cross-donor ASM claim

**The test is correctly calibrated — the heterogeneity is biological.** Empirical null (split one
haplotype's reads in half, run the identical test), 199 donors: **λ_null = 0.61–0.71, zero donors
above 1.1**, i.e. uniformly slightly conservative.

**But donors differ enormously in real haplotype divergence.** In the real comparison, chr20
λ ranges 0.81–13.1 (median 1.54), so raw ASM yield spans 0.02%–11% per donor.

- λ tracks **domain depth** (p < 1e-4) and **clonality** (r = 0.40), not chemistry (p = 0.67) or
  ancestry (p > 0.2 with depth in the model). The apparent 4x ancestry difference in raw ASM
  yield is mediated by donor state.
- λ is uniform across chromosomes within a donor, so it is not aneuploidy or LOH.
- **Excluding PMD-overlapping calls does not fix it**: λ outside domains is nearly as high
  (NA20762 10.2 vs 14.6 inside; NA19338 higher outside).
- Mechanism: within-haplotype read spread correlates *negatively* with λ (r = −0.705). In clonal
  or domain-deep lines the two alleles' states are not averaged across cells, so real allelic
  differences appear genome-wide — drift, not locus-specific regulation.

**Handling (do not filter on clonality — it is female-only and insufficient; 44 inflated donors
remain after dropping all 35 clonal-like ones):**

| tier (chr20 λ) | donors | handling |
|---|---|---|
| < 1.2 | 70 | use as-is |
| 1.2–3 | 106 | genomic control or empirical null |
| > 3 | 25 | exclude from genome-wide discovery; keep for targeted tests |
| > 8 | ~3 (incl. NA20762) | exclude |

Genomic control preserves imprinted-DMR detection up to λ ≈ 6–8 (HG00097 113→113, NA18508
107→101, HG01150 94→76) and over-corrects only at λ = 12.9 (NA20762 91→0).

**Consequence for penetrance:** with 66% of union loci private to one donor and per-donor yield
scaling with λ, penetrance is only interpretable after calibration. Recipe: discover on λ < 1.2
donors, re-test candidates in all donors with genomic control, define penetrance over *tested and
calibrated* donors, and anchor against imprinted DMRs and HPRC2 mQTLs.

---

## Where the data lives

| what | path |
|---|---|
| Donor × CpG matrices (union CpGs; `n_meth`, `n_total`; all + het-filtered) | `results/asm/cpg_matrix/<chrom>.{all,hetfilt}.npz` — 30.9M CpGs × 402 haplotypes |
| ASM calls / genome-wide merge | `results/asm/calls/`, `results/asm/genome/` |
| ASM empirical null | `results/asm/null/`, `results/qc/data/asm_null_vs_real_chr20.tsv` |
| Genome-wide PMD calls (native) | `data/pmds/<s>/<s>_hap<h>.pmd.bed` |
| 10 kb bin matrix, domain frequency, PCA | `results/meth_bins/` |
| RT / LAD annotation | `results/meth_bins/annotations/rt_lad_10kb.tsv.gz` |
| Solo-WCGW | `results/meth_bins/solo_wcgw/` |
| Per-CpG hg38 bigWigs | `data/pmds/<s>/<s>_hap<h>.hg38.{meth,depth}.bw` |
| UCSC hubs | `/u/project/cluo/PUBLIC_SHARED/ucsc/asm_lr_hprc2{,_gradient}` |
| Figures | `figures/{qc,pmds,asm,genetics}/` |
