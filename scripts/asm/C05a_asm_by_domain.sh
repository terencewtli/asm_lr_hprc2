#!/bin/bash
#$ -N C05a_asm_by_domain
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 4
#$ -hold_jid C04a_merge_asm_genome
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C05a_asm_by_domain.$JOB_ID
#$ -j y

# Is ASM concentrated in PMDs, and does it look genetic (recurrent across donors, mQTL-linked)
# or stochastic (donor-private)? Held on the genome-wide ASM merge.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/asm"
echo "$(date): C05a"
time python3 C05a_asm_by_domain.py
echo "$(date): done"
