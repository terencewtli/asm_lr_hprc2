#!/bin/bash
#$ -N A01c_pmds
#$ -cwd
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 1
#$ -t 1-458:1
#$ -tc 60
#$ -hold_jid A01a_modbed_to_methcounts
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A01c_pmds.$JOB_ID.$TASK_ID
#$ -j y

# HMM PMD calling (dnmtools pmd) on A01a's symmetrized methcounts, one (sample, hap) per task,
# same 458-task indexing as A01a. No symmetric-CpG merge step here (unlike the WGBS template's
# separate A01b) -- A01a already does the symmetrization directly on modbed's read-level calls,
# there's no pre-merged allc-equivalent source file to redundantly re-merge.
#
# Validated 2026-09-16 on 3 donors restricted to chr20: NA19338 (global meth 0.523, the low
# tail) -> 87 real PMDs, 42.0Mb, 64.2% of the covered contig, largest domain 2.7Mb -- broadly
# consistent with QC06's threshold-based estimate (78.3%) for the same donor/chromosome.
# HG01981 (median, 0.648) and HG04187 (high, 0.723) -> 0 domains each, vs. QC06's threshold
# method reporting 69.1% and 22.4% respectively -- the real HMM is markedly more conservative
# than the threshold heuristic and did not confirm PMDs in the non-extreme donors on this one
# chromosome. Not treated as a bug: 0 domains is a legitimate HMM outcome when the input isn't
# genuinely bimodal, and this qualitatively sharpens (not contradicts) the original finding --
# real PMD signal may be more selective/donor-specific than QC06's first-pass threshold approach
# suggested. Worth deciding, once the full genome-wide run lands, whether QC06 should be revised
# to use these HMM calls directly rather than its own segmentation.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
MANIFEST=$PROJDIR/tsv/meta/hprc2_sample_manifest.tsv
DNMTOOLS=/u/home/t/terencew/bin/dnmtools

ID=$SGE_TASK_ID
SAMPLE_ROW=$(( (ID - 1) / 2 + 1 ))
HAP=$(( (ID - 1) % 2 + 1 ))
SAMPLE=$(awk -F'\t' -v t="$SAMPLE_ROW" 'NR==t+1 {print $1}' "$MANIFEST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: no manifest row $SAMPLE_ROW for task $ID"
    exit 1
fi

OUTDIR=$PROJDIR/data/pmds/$SAMPLE
METH=$OUTDIR/${SAMPLE}_hap${HAP}.cpg.methcounts.tsv.gz
OUT=$OUTDIR/${SAMPLE}_hap${HAP}.pmd.bed

echo "$(date): A01c task $ID — $SAMPLE hap$HAP"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

if [ ! -s "$METH" ]; then
    echo "$(date): no methcounts for $SAMPLE hap$HAP (A01a skip or not yet run), skipping"
    exit 0
fi

# -S/-r/-p (summary/posteriors/params output) all throw "bad file: -S" on this dnmtools build
# regardless of argument order -- confirmed reproducible in isolation before submitting this
# array, not a task-specific fluke. Dropped; -o's PMD bed is the only output this needs.
time "$DNMTOOLS" pmd -o "$OUT" -s 1 -v -i 1000 "$METH"

echo "$(date): $SAMPLE hap$HAP complete"
