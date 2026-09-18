#!/bin/bash
#$ -N C07d_replication_summary
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 10
#$ -hold_jid C07c_candidate_genotypes
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C07d_replication_summary.$JOB_ID
#$ -j y

# ASM replication step 4: penetrance on replication-tier donors, replication classes, and the
# replicating-vs-not feature comparison. Held on the whole C07c array.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/asm"
echo "$(date): C07d"
time python3 C07d_replication_summary.py
echo "$(date): done"
