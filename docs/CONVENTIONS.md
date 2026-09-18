# Conventions — how this project is operated

Read with `PROGRESS.md` (what has run) and `RESULTS.md` (what we know). This file is the
operating manual: the things a new session has to know before touching anything.

## Environment

- Python: the **`allcools`** conda env (`/u/home/t/terencew/project-cluo/miniconda3/envs/allcools`).
  It is **Python 3.7** — no `math.comb`, no `Path.unlink(missing_ok=)`, no walrus-free-ness
  assumptions from 3.8+. `claude_env` is Node.js, not Python.
- R: kernel `ir410`. In a qsub job you must
  `source /u/local/Modules/default/init/modules.sh; module load R/4.1.0` first, or the IRkernel
  dies on missing MKL libraries.
- Tools: `dnmtools`, `bedtools`, `bigWig*` in `~/bin`; UCSC kent utils in
  `/u/project/cluo/terencew/programs/kentsrc/utils`; `liftOver` is there but is too slow for
  per-CpG work — use `scripts/meth_bins/B01a_hap_bins.map_points` instead.

## Script naming

`<Letter><NN><letter>_<name>.{py,sh}` — the letter is the pipeline, the number is the task group,
the trailing letter is the step:

| prefix | pipeline |
|---|---|
| `A01*`, `A02*` | modbed → methcounts → PMD calls → bigWigs → metagene (`scripts/call_pmds/`, `scripts/bigwig/`) |
| `B01*`–`B06*` | 10kb bin matrix, domain frequency, QC, variant density, solo-WCGW, RT/LAD (`scripts/meth_bins/`) |
| `C01*`–`C07*` | ASM: regions → per-donor calls → CpG matrices → genome merge → domain join → empirical null → replication (C07a candidates, b re-test, c genotype context, d classes) (`scripts/asm/`) |
| `G0*`, `H0*`, `P0*` | assembly-vs-hg38 variants, het sites, ancestry PCA (`scripts/harmonize_hg38/`, `scripts/qsub/`) |
| `M0*`, `Q0*`, `R0*` | one-off per-donor metrics: global mCG, XIST skew, molecule QC, RNA markers |
| `U01*` | UCSC track hubs (`scripts/ucsc/`) |

Notebooks: `notebooks/<analysis type>/<NN><letter>_<name>.ipynb`; figure notebooks in
`notebooks/final_figures/figure_<n|sN>/{python,R}` — python exports CSVs to `csv/figure_*`, R
reads those and writes PDFs to `pdf/figure_*`. Figures from analysis notebooks go to
`figures/<category>/`.

## Array-job pattern (follow it, it has bitten us)

Every array script: `#$ -t 1-N`, a `# ID=1` commented line above `ID=$SGE_TASK_ID` so the script
is testable by hand, skip-if-exists on its own output, an explicit skip (exit 0) for donors with
missing inputs, `time` before the main command, logs to `logs/<jobname>.$JOB_ID.$TASK_ID`.

**The skip-if-exists check lives in the shell script, not the Python.** If you change the Python
output path, change the shell check too — otherwise every task silently skips and the job
"succeeds" having done nothing (this cost a full G02 rerun).

Task lists: `txt/samples/h01_sample_chrom_tasks.txt` (4,444 sample×chrom rows) and the manifest
`tsv/meta/hprc2_sample_manifest.tsv` (229 rows; only 202 have assemblies — index by manifest row,
never assume the first 202 rows are the usable ones).

## Data conventions

- Coordinates: PMDs are called in each haplotype's **own assembly coordinates**; per-CpG values
  are mapped to **hg38** for anything cross-donor. Autosomes only.
- A CpG is only counted if the assembly reads `CG` there; strands are merged onto the `+`
  coordinate; strand-flipped chain blocks need the −1 correction.
- Every cross-donor model must include **ONT chemistry**; every ASM comparison must account for
  **per-donor inflation (λ)**. See RESULTS §1 and §7.

## Before trusting a directory

Check `PROGRESS.md`'s "Known-stale / do-not-use" list first. Stale directories also carry a
`DEPRECATED.md` file explaining what replaced them.
