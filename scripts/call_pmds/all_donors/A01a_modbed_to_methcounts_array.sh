#!/bin/bash
#$ -N A01a_modbed_to_methcounts
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 1
#$ -t 1-458:1
#$ -tc 60
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A01a_modbed_to_methcounts.$JOB_ID.$TASK_ID
#$ -j y

# Whole-genome (every contig in the modbed, no chr restriction) modbed -> symmetrized methcounts
# conversion, one (sample, hap) per task -- 229 manifest rows x 2 haps = 458 tasks. Replaces the WGBS
# template's A01a_allc_to_methcounts.sh; see A01a_modbed_to_methcounts.py's docstring for why
# this needed real per-position aggregation + CpG symmetrization, not a column rename, and for
# the empirical check that resolved the strand/symmetrization question this session's earlier
# validation run got wrong on the first pass (dnmtools pmd requires strand-collapsed input;
# confirmed via its own usage text, not assumed).
#
# Validated 2026-09-16 on 3 donors (HG00097, NA19338, HG04187) restricted to chr20 before this
# full run was submitted -- see A01c's array script for the PMD-calling step this holds for.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
MANIFEST=$PROJDIR/tsv/meta/hprc2_sample_manifest.tsv

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID

# 458 tasks = 229 manifest rows x 2 haps, hap-major within each sample (task 1,2 = sample1
# hap1,hap2; task 3,4 = sample2 hap1,hap2; ...)
SAMPLE_ROW=$(( (ID - 1) / 2 + 1 ))
HAP=$(( (ID - 1) % 2 + 1 ))
SAMPLE=$(awk -F'\t' -v t="$SAMPLE_ROW" 'NR==t+1 {print $1}' "$MANIFEST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: no manifest row $SAMPLE_ROW for task $ID"
    exit 1
fi

OUTDIR=$PROJDIR/data/pmds/$SAMPLE
mkdir -p "$OUTDIR"
OUT=$OUTDIR/${SAMPLE}_hap${HAP}.cpg.methcounts.tsv.gz

echo "$(date): A01a task $ID — $SAMPLE hap$HAP"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

if [ ! -s "$PROJDIR/data/modbed/${SAMPLE}_hap${HAP}.modbed.gz" ]; then
    echo "$(date): no modbed for $SAMPLE hap$HAP, skipping (not a task failure)"
    exit 0
fi

# The reference-CpG filter needs the haplotype assembly; 27/229 manifest donors have none
# (unresolved assembly naming), so they are skipped rather than written unfiltered.
if [ ! -s "$PROJDIR/data/assemblies/${SAMPLE}_hap${HAP}.fa.gz" ]; then
    echo "$(date): no assembly for $SAMPLE hap$HAP, skipping (not a task failure)"
    exit 0
fi

cd "$PROJDIR/scripts/call_pmds/all_donors"
time python3 A01a_modbed_to_methcounts.py "$SAMPLE" "$HAP" "$OUT"

echo "$(date): $SAMPLE hap$HAP complete"
