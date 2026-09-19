#!/bin/bash
#$ -N C08c_lead_permutation
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 10
#$ -t 1-22:1
#$ -tc 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C08c_lead_permutation.$JOB_ID.$TASK_ID
#$ -j y

# Permutation null for the C07c lead-variant scan: one task per autosome, 1000 permutations of
# which calibrated donors carry the ASM call, holding the variant set / het matrix / K fixed.
# Corrects the uncorrected P_LEAD = 1e-4 threshold in C07d.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2

# ID=20
ID=$SGE_TASK_ID

CHROM=chr$ID
OUT=$PROJDIR/results/asm/model/lead_perm/$CHROM.lead_perm.tsv.gz
if [ -s "$OUT" ]; then
    echo "$(date): exists, skipping $OUT"
    exit 0
fi
echo "$(date): C08c $CHROM"
cd "$PROJDIR/scripts/asm"
time python3 C08c_lead_permutation.py "$CHROM" 1000
echo "$(date): done"
