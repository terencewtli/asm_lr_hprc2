#!/bin/bash
#$ -N M03_xist_promoter_skew
#$ -cwd
#$ -l h_data=3G,h_rt=2:00:00
#$ -pe shared 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/M03_xist_promoter_skew.$JOB_ID
#$ -j y

# XCI-skew / clonality metric: per-haplotype methylation of the XIST promoter CpG island in every
# female donor (QC05's gene-body window was too diluted; see script docstring).

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/qsub"
echo "$(date): M03"
time python3 M03_xist_promoter_skew.py 10
echo "$(date): done"
