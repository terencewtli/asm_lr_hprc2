#!/bin/bash
#$ -N B01a_hap_bins
#$ -cwd
#$ -l h_data=4G,h_rt=2:00:00
#$ -pe shared 4
#$ -t 1-458:1
#$ -tc 60
#$ -hold_jid A01c_pmds
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/B01a_hap_bins.$JOB_ID.$TASK_ID
#$ -j y

# Per-(sample, hap) genome-wide hg38 10kb methylation bins + hg38 PMD intervals from A01a's
# CpG methcounts and A01c's genome-wide PMD calls. Same 229x2=458 indexing as A01a/A01c; held on
# A01c. Aggregated by B01b_bin_matrix_pca.py. 4 slots x 4G = memory headroom for ~32M CpGs in pandas.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
MANIFEST=$PROJDIR/tsv/meta/hprc2_sample_manifest.tsv

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID
SAMPLE_ROW=$(( (ID - 1) / 2 + 1 ))
HAP=$(( (ID - 1) % 2 + 1 ))
SAMPLE=$(awk -F'\t' -v t="$SAMPLE_ROW" 'NR==t+1 {print $1}' "$MANIFEST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: no manifest row $SAMPLE_ROW for task $ID"
    exit 1
fi

OUT=$PROJDIR/results/meth_bins/per_hap/${SAMPLE}_hap${HAP}.summary.tsv
echo "$(date): B01a task $ID — $SAMPLE hap$HAP"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

for f in "data/pmds/$SAMPLE/${SAMPLE}_hap${HAP}.cpg.methcounts.tsv.gz" \
         "data/pmds/$SAMPLE/${SAMPLE}_hap${HAP}.pmd.bed" \
         "data/chains/${SAMPLE}_hap${HAP}_vs_GRCh38.chain.gz"; do
    if [ ! -e "$PROJDIR/$f" ]; then
        echo "$(date): missing $f, skipping (not a task failure)"
        exit 0
    fi
done

cd "$PROJDIR/scripts/meth_bins"
time python3 B01a_hap_bins.py "$SAMPLE" "$HAP"

echo "$(date): $SAMPLE hap$HAP complete"
