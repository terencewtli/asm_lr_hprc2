#!/bin/bash
#$ -N C03a_cpg_matrix
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 10
#$ -t 1-22:1
#$ -hold_jid C02a_call_asm,A02a_methcounts_to_bigwig
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C03a_cpg_matrix.$JOB_ID.$TASK_ID
#$ -j y

# Per-chromosome hg38 donor-haplotype x CpG matrices over the union of CpGs, both sources
# (hetfilt from C02a, all from the A02a bigWigs). Task ID = chromosome number.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2

# ID=20  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
CHROM=chr$ID

cd "$PROJDIR/scripts/asm"
for SRC in hetfilt all; do
    OUT=$PROJDIR/results/asm/cpg_matrix/${CHROM}.${SRC}.npz
    if [ -s "$OUT" ]; then
        echo "$(date): $OUT exists, skipping"
        continue
    fi
    echo "$(date): C03a $CHROM $SRC"
    time python3 C03a_cpg_matrix.py "$CHROM" "$SRC"
done
echo "$(date): $CHROM complete"
