#!/bin/bash
#$ -N C07c_candidate_genotypes
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 10
#$ -t 1-22:1
#$ -hold_jid C07b_retest_candidates
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C07c_candidate_genotypes.$JOB_ID.$TASK_ID
#$ -j y

# ASM replication step 3: one task per autosome. Genotype context per (donor, candidate) and the
# per-candidate lead-variant association. Held on the whole C07b array (needs every donor's calls).

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2

# ID=20
ID=$SGE_TASK_ID

CHROM=chr$ID
OUT=$PROJDIR/results/asm/replication/genotype/$CHROM.region.tsv.gz
if [ -s "$OUT" ]; then
    echo "$(date): exists, skipping $OUT"
    exit 0
fi
echo "$(date): C07c $CHROM"
cd "$PROJDIR/scripts/asm"
time python3 C07c_candidate_genotypes.py "$CHROM"
echo "$(date): done"
