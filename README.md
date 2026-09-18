# asm_lr_hprc2 — population-scale allele-specific methylation called from individuals in HPRC2

## Original project conception

Successor to [`asm_lr`](https://github.com/terencewtli/asm_lr), re-platformed onto the HPRC
Release 2 (HPRC2) epigenome resource (`epigenome.humanpangenome.org`, S3 bucket
`hprc-epigenome`). Check documentation in `docs/20260908_project_transition.md`.

## Core question

| Superpop | n |
|---|---|
| AFR | 63 |
| AMR | 43 |
| EAS | 49 |
| EUR | 30 |
| SAS | 36 |
| unmapped | 8 |

Allele-specific methylation (ASM) has been studied in both short-read/long-read contexts, but the mechanisms of ASM penetrance
have not been elucidated at a population genetics level. Penetrance in this context refers to the proportion of a particular haplotype that
exhibits ASM. Additionally, mechanisms such as allele frequencies, CpG architecture, and
LD can themselves lead to an apparent lack of ASM (statistical vs. biological impact on ASM calling) between populations - these design features
have previously not been appreciated across multiple populations when considering ASM.
Leveraging ancestrally diverse haplotypes, our goal is to uncover the sequence determinants of ASM to better understand this form of regulation.

## Pipeline

HPRC2 Epigenome browser contains donor hap1/hap2 modbed files aligned to each donor's haploid assembly:

https://epigenome.humanpangenome.org/?tab=sample

Using HPRC1 data, we have developed a beta-binomial, coverage-aware within-sample ASM caller
to derive ASM calls from ONT data. https://github.com/terencewtli/ont_asm_caller

Full schematic and details of pipeline are pending.

## Data & infrastructure

- `data/hprc2_sample_manifest.tsv` — all 229 usable samples: population/superpopulation, ONT
  modbed file sizes (coverage proxy), other available data types, resolved assembly download
  URLs where available, overlap flags vs. `asm_lr`'s 18- and 30-donor cohorts.
- `docs/data_sources.md` — how every download link and format claim was verified.

## Related

- [`hprc2_misc`](https://github.com/terencewtli/hprc2_misc) — sister repo: HPRC2 RNA / Hi-C /
  Fiber-seq for these donors (orthogonal PMD validation), de novo meQTL mapping, and other
  analyses outside the ASM-penetrance question.
