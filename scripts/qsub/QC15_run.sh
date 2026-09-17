#!/bin/bash
#$ -N QC15_run
#$ -cwd
#$ -l h_data=6G,h_rt=4:00:00
#$ -pe shared 6
#$ -hold_jid B04a_variant_density,B05a_solo_wcgw,M03_xist_promoter_skew
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/QC15_run.$JOB_ID
#$ -j y

# Execute QC15 once RT/LAD, variant density, solo-WCGW and XIST-skew inputs all exist.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/notebooks/pmds"
echo "$(date): QC15"
time jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=7200 \
    --ExecutePreprocessor.kernel_name=allcools 04a_mechanism_clock_instability.ipynb
echo "$(date): done"
