#!/bin/bash
#$ -N B04a_variant_density
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 2
#$ -t 1-22:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/B04a_variant_density.$JOB_ID.$TASK_ID
#$ -j y

# Per-10kb variant density (SNV / indel / >=50bp SV) from G01's per-chrom assembly-vs-hg38
# diploid VCFs, one chromosome per task, summed across ~200 donors.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
# ID=20  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
CHROM=chr$ID

OUT=$PROJDIR/results/meth_bins/variant_density/${CHROM}.bins.tsv.gz
if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

cd "$PROJDIR/scripts/meth_bins"
echo "$(date): B04a $CHROM"
time python3 B04a_variant_density.py "$CHROM"
echo "$(date): done"
