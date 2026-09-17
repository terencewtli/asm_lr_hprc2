#!/bin/bash
#$ -N P01_1000g_full_pca
#$ -cwd
#$ -l h_data=8G,h_rt=12:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/P01_1000g_full_pca.$JOB_ID
#$ -j y

# Classic 1000G PCA on the FULL 3,202-sample panel, so HPRC2 donors (themselves 1000G samples)
# can be shown against the whole reference cohort rather than against each other only.
# Reuses the shared script; output lands in the shared reference tree, not this project.

set -euo pipefail
export PATH="/u/local/apps/bcftools/1.11/gcc-4.8.5/bin:$PATH"
REF=/u/project/cluo/terencew/claude/reference/1000G
bash $REF/scripts/run_1000g_pca.sh $REF/tsv/meta/g1k_3202_sample_ids.txt $REF/pca/g1k_full
echo "$(date): done"
