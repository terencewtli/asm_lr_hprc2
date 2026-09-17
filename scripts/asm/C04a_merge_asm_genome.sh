#!/bin/bash
#$ -N C04a_merge_asm_genome
#$ -cwd
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 10
#$ -hold_jid C02a_call_asm
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C04a_merge_asm_genome.$JOB_ID
#$ -j y

# Genome-wide ASM merge (per-donor genome-wide BH, region penetrance by chemistry, Zink
# imprinted-DMR control by chemistry). Held on the whole C02a array.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/asm"
echo "$(date): C04a"
time python3 C04a_merge_asm_genome.py
echo "$(date): done"
