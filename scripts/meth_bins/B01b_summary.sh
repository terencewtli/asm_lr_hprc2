#!/bin/bash
#$ -N B01b_summary
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 10
#$ -hold_jid B01a_hap_bins
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/B01b_summary.$JOB_ID
#$ -j y

# Genome-wide aggregation once every B01a task is done: sample x 10kb bin matrix + PCA +
# covariate R^2 + continuum metrics + domain frequency (B01b), then LCL vs fibroblast PMD
# boundaries/genes (B02a). Both write under results/meth_bins/ and are safe to rerun.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/meth_bins"

echo "$(date): B01b bin matrix / PCA"
time python3 B01b_bin_matrix_pca.py
echo "$(date): B02a boundary genes"
time python3 B02a_pmd_boundary_genes.py
echo "$(date): done"
