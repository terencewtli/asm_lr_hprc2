#!/bin/bash
#$ -N C10a_epiallele_nme
#$ -cwd
#$ -l h_data=6G,h_rt=8:00:00
#$ -pe shared 10
#$ -t 1-69:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/C10a_epiallele_nme.$JOB_ID.$TASK_ID
#$ -j y

# Epiallele entropy (CPEL NME) at the 10,756 replicating ASM loci, discovery-tier donors only.
# Scoped deliberately: NME + MML only (no PDM), no locus calling. Pre-registered question is
# whether epiallele entropy predicts penetrance (RESULTS §13). Empirical null is built in.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
LIST=$PROJDIR/txt/samples/c10a_discovery_donors.txt

# ID=1
ID=$SGE_TASK_ID

SAMPLE=$(sed -n "${ID}p" "$LIST")
OUT=$PROJDIR/results/asm/model/epiallele/$SAMPLE.nme.tsv.gz
if [ -s "$OUT" ]; then
    echo "$(date): exists, skipping $OUT"
    exit 0
fi
echo "$(date): C10a $SAMPLE"
cd "$PROJDIR/scripts/asm"
time python3 C10a_epiallele_nme.py "$SAMPLE"
echo "$(date): done"
