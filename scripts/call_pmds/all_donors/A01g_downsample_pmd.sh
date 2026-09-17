#!/bin/bash
#$ -N A01g_downsample_pmd
#$ -cwd
#$ -l h_data=4G,h_rt=6:00:00
#$ -pe shared 4
#$ -t 1-4:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A01g_downsample_pmd.$JOB_ID.$TASK_ID
#$ -j y

# Coverage-sufficiency test for PMD calling: thin a haplotype's methcounts to 20/15/10/5x and
# re-call. Four tasks = two donors spanning the domain-depth range, hap1 and hap2.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail
PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
case $ID in
  1) SAMPLE=NA19338; HAP=1 ;;   # deep domains
  2) SAMPLE=NA19338; HAP=2 ;;
  3) SAMPLE=HG00097; HAP=1 ;;   # typical
  4) SAMPLE=HG00097; HAP=2 ;;
esac

cd "$PROJDIR/scripts/call_pmds/all_donors"
echo "$(date): A01g $SAMPLE hap$HAP"
time python3 A01g_downsample_pmd.py "$SAMPLE" "$HAP"
echo "$(date): done"
