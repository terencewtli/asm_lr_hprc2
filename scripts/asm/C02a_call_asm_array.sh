#!/bin/bash
#$ -N C02a_call_asm
#$ -cwd
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-4444:1
#$ -tc 100
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C02a_call_asm.$JOB_ID.$TASK_ID
#$ -j y

# Read-level ASM calls + het-filtered hg38 per-CpG counts for one (sample, chrom); one task per
# row of H01's task list (202 donors x chr1-22). Replaces P03_build_donor_chrom_matrix_array.sh.
# Test (2026-09-17): HG00097 chr20, 2.3 min, 2.8 GB RSS, 49,842 regions tested.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
TASKS=$PROJDIR/txt/samples/h01_sample_chrom_tasks.txt

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
SAMPLE=$(awk -F'\t' -v t="$ID" 'NR==t {print $1}' "$TASKS")
CHROM=$(awk -F'\t' -v t="$ID" 'NR==t {print $2}' "$TASKS")

OUT=$PROJDIR/results/asm/calls/${SAMPLE}_${CHROM}.asm.tsv.gz
echo "$(date): C02a task $ID — $SAMPLE $CHROM"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

for f in "data/chains/${SAMPLE}_hap1_vs_GRCh38.chain.gz" "data/chains/${SAMPLE}_hap2_vs_GRCh38.chain.gz" \
         "data/het_snps/tmp_${CHROM}/${SAMPLE}/het_in_hap1.bed" "data/hmmflagger/${SAMPLE}_HMMFlagger.ONT.bed.gz"; do
    if [ ! -s "$PROJDIR/$f" ]; then
        echo "$(date): missing $f, skipping (not a task failure)"
        exit 0
    fi
done

cd "$PROJDIR/scripts/asm"
time python3 C02a_call_asm.py "$SAMPLE" "$CHROM"

echo "$(date): $SAMPLE $CHROM complete"
