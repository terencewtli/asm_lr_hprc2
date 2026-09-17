#!/bin/bash
#$ -N B05a_solo_wcgw
#$ -cwd
#$ -l h_data=3G,h_rt=12:00:00
#$ -pe shared 10
#$ -hold_jid A02a_methcounts_to_bigwig
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/B05a_solo_wcgw.$JOB_ID
#$ -j y

# Solo-WCGW (Zhou 2018 mitotic clock) methylation per haplotype, by LCL domain class.
# ~5 min per haplotype; skips finished haplotypes. Sites annotated once (sites.npz).

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
cd "$PROJDIR/scripts/meth_bins"
echo "$(date): B05a"
time python3 B05a_solo_wcgw.py run 10
echo "$(date): done"
