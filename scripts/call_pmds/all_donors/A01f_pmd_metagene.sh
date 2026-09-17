#!/bin/bash
#$ -N A01f_pmd_metagene
#$ -cwd
#$ -l h_data=2G,h_rt=8:00:00
#$ -pe shared 10
#$ -hold_jid A02a_methcounts_to_bigwig
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A01f_pmd_metagene.$JOB_ID
#$ -j y

# PMD metagene profiles for every haplotype bigWig (A02a) over the A01e reference sets.
# ~1 min per haplotype single-threaded (HG00097 hap1 test, 2026-09-17); skips finished haplotypes.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/call_pmds/all_donors"
echo "$(date): A01f"
time python3 A01f_pmd_metagene.py 10
echo "$(date): done"
