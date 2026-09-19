#!/bin/bash
#$ -N C09a_haplotype_background
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00
#$ -pe shared 10
#$ -t 1-22:1
#$ -tc 10
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C09a_haplotype_background.$JOB_ID.$TASK_ID
#$ -j y

# Does ASM penetrance depend on WHICH local haplotype carries the lead allele (RESULTS §13)?
# One task per autosome, 1000 permutations of which lead-het carriers show ASM.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2

# ID=21
ID=$SGE_TASK_ID

CHROM=chr$ID
OUT=$PROJDIR/results/asm/model/hap_background/$CHROM.hapbg.tsv.gz
if [ -s "$OUT" ]; then
    echo "$(date): exists, skipping $OUT"
    exit 0
fi
echo "$(date): C09a $CHROM"
cd "$PROJDIR/scripts/asm"
time python3 C09a_haplotype_background.py "$CHROM" 1000
echo "$(date): done"
