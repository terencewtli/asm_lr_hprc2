#!/bin/bash
#$ -N U01c_rebuild_hubs_and_qc12
#$ -cwd
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 10
#$ -hold_jid A02a_methcounts_to_bigwig,A01f_pmd_metagene
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/U01c_rebuild_hubs_and_qc12.$JOB_ID
#$ -j y

# After all per-CpG bigWigs (A02a) and metagene profiles (A01f) exist: rebuild both UCSC hubs
# (adds the per-CpG tracks that were missing on the first build) and re-execute QC12 so its
# metagene section covers every haplotype.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/ucsc"
echo "$(date): U01a main hub"
time python3 U01a_build_track_hub.py
echo "$(date): U01b gradient hub"
time python3 U01b_build_gradient_hub.py
echo "$(date): QC12 notebook"
cd "$PROJDIR/notebooks/pmds"
time jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=7200 \
    --ExecutePreprocessor.kernel_name=allcools 03a_pmd_size_overlap_and_metagene.ipynb
echo "$(date): done"
