#!/bin/bash
#$ -N A03a_chr20_pmd_survey
#$ -cwd
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-458:1
#$ -tc 100
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A03a_chr20_pmd_survey.$JOB_ID.$TASK_ID
#$ -j y

# chr20-only PMD survey per (sample, hap): CpG-filtered methcounts -> dnmtools pmd -> hg38-lifted
# per-CpG table. Same 229x2=458 indexing as A01a/A01c. Independent of the genome-wide A01a run
# (reads modbed directly), so no hold. ~2 min/task (NA19338 hap1 test, 2026-09-16).
# Aggregated by A03b_chr20_pmd_survey_summary.py.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
MANIFEST=$PROJDIR/tsv/meta/hprc2_sample_manifest.tsv
TMPROOT=/u/project/cluo_scratch/terencew/claude/asm_lr_hprc2/chr20_survey

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
SAMPLE_ROW=$(( (ID - 1) / 2 + 1 ))
HAP=$(( (ID - 1) % 2 + 1 ))
SAMPLE=$(awk -F'\t' -v t="$SAMPLE_ROW" 'NR==t+1 {print $1}' "$MANIFEST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: no manifest row $SAMPLE_ROW for task $ID"
    exit 1
fi

OUT=$PROJDIR/results/qc/data/chr20_pmd_survey/per_hap/${SAMPLE}_hap${HAP}.summary.tsv
echo "$(date): A03a task $ID — $SAMPLE hap$HAP"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

for f in "data/modbed/${SAMPLE}_hap${HAP}.modbed.gz" "data/assemblies/${SAMPLE}_hap${HAP}.fa.gz" \
         "data/chains/${SAMPLE}_hap${HAP}_vs_GRCh38.chain.gz"; do
    if [ ! -s "$PROJDIR/$f" ]; then
        echo "$(date): missing $f, skipping (not a task failure)"
        exit 0
    fi
done

cd "$PROJDIR/scripts/call_pmds/chr20_survey"
time python3 A03a_chr20_pmd_survey.py "$SAMPLE" "$HAP" "$TMPROOT/${SAMPLE}_hap${HAP}"
rmdir "$TMPROOT/${SAMPLE}_hap${HAP}" 2>/dev/null || true

echo "$(date): $SAMPLE hap$HAP complete"
