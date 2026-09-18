#!/bin/bash
#$ -N S1_figures_run
#$ -cwd
#$ -l h_data=6G,h_rt=4:00:00
#$ -pe shared 4
#$ -hold_jid Q01a_molecule_qc,P01_1000g_full_pca
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/S1_figures_run.$JOB_ID
#$ -j y

# Supplementary Figure S1 (cohort QC): python notebook assembles the per-donor/per-haplotype QC
# table and exports CSVs, then the R notebook draws the panels. Held on the molecule-level QC
# array and the full 1000G PCA (panel I falls back to the 221-donor PCA if that isn't ready).

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2

echo "$(date): S1 python"
cd "$PROJDIR/notebooks/final_figures/figure_s1/python"
time jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=7200 \
    --ExecutePreprocessor.kernel_name=allcools figure_s1_qc_python.ipynb

# R 4.1.0 links against MKL; without the module its shared libs are missing and the IRkernel
# dies with "Kernel died before replying to kernel_info" (cost one run, 2026-09-17).
source /u/local/Modules/default/init/modules.sh
module load R/4.1.0

echo "$(date): S1 R"
cd "$PROJDIR/notebooks/final_figures/figure_s1/R"
time jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=7200 \
    --ExecutePreprocessor.kernel_name=ir410 figure_s1_R.ipynb
echo "$(date): done"
