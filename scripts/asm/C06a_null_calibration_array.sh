#!/bin/bash
#$ -N C06a_null_calibration
#$ -cwd
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-202:1
#$ -tc 100
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C06a_null_calibration.$JOB_ID.$TASK_ID
#$ -j y

# Empirical null for the ASM test: split one haplotype's reads in half and run the same test.
# chr20 for every donor (one task per donor). See C06a_null_calibration.py docstring.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
TASKS=$PROJDIR/txt/samples/h01_sample_chrom_tasks.txt
CHROM=chr20

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
SAMPLE=$(awk -F'\t' -v c="$CHROM" '$2==c {print $1}' "$TASKS" | sed -n "${ID}p")
if [ -z "$SAMPLE" ]; then echo "no sample for task $ID"; exit 0; fi

OUT=$PROJDIR/results/asm/null/${SAMPLE}_${CHROM}.null_summary.tsv
if [ -s "$OUT" ]; then echo "$(date): $OUT exists, skipping"; exit 0; fi
if [ ! -s "$PROJDIR/data/chains/${SAMPLE}_hap1_vs_GRCh38.chain.gz" ]; then
    echo "$(date): no chain for $SAMPLE, skipping (not a task failure)"; exit 0
fi

cd "$PROJDIR/scripts/asm"
echo "$(date): C06a $SAMPLE $CHROM"
time python3 C06a_null_calibration.py "$SAMPLE" "$CHROM"
echo "$(date): done"
