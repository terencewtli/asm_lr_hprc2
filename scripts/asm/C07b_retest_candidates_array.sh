#!/bin/bash
#$ -N C07b_retest_candidates
#$ -cwd
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 10
#$ -t 1-201:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C07b_retest_candidates.$JOB_ID.$TASK_ID
#$ -j y

# ASM replication step 2: one task per donor in results/asm/replication/donor_tiers.tsv (C07a),
# re-tests the candidate loci with per-donor genomic control. Skip-if-exists on the output.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
TIERS=$PROJDIR/results/asm/replication/donor_tiers.tsv

# ID=1
ID=$SGE_TASK_ID

SAMPLE=$(awk -F'\t' -v t="$ID" 'NR==t+1 {print $1}' "$TIERS")
if [ -z "$SAMPLE" ]; then
    echo "no donor for task $ID"
    exit 0
fi
OUT=$PROJDIR/results/asm/replication/calls/$SAMPLE.tsv.gz
if [ -s "$OUT" ]; then
    echo "$(date): exists, skipping $OUT"
    exit 0
fi
echo "$(date): C07b $SAMPLE"
cd "$PROJDIR/scripts/asm"
time python3 C07b_retest_candidates.py "$SAMPLE"
echo "$(date): done"
