#!/bin/bash
#$ -N QC16_run
#$ -cwd
#$ -l h_data=6G,h_rt=4:00:00
#$ -pe shared 10
#$ -hold_jid C06a_null_calibration
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/QC16_run.$JOB_ID
#$ -j y

# ASM calibration/QC dashboard; held until the empirical-null array (C06a) finishes.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/notebooks/asm"
echo "$(date): QC16"
time jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=7200 \
    --ExecutePreprocessor.kernel_name=allcools 01a_asm_calibration_qc.ipynb
echo "$(date): done"
