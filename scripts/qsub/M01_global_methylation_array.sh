#!/bin/bash
#$ -N M01_global_methylation_array
#$ -cwd
#$ -l h_data=2G,h_rt=0:30:00
#$ -pe shared 1
#$ -t 1-229:1
#$ -tc 60
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/M01_global_methylation_array.$JOB_ID.$TASK_ID
#$ -j y

# One array task per manifest row. Computes REAL genome-wide (not subsampled) global methylation
# for both haplotypes of that donor -- full zcat+awk scan of modbed col7 (meth offsets) / col8
# (unmeth offsets) per docs/data_sources.md SS5, same counting logic as the n=10 validation in
# notebooks/qc/QC04_global_methylation_covariates.ipynb, just moved from a 1-core interactive
# shell (~90s/file serial) to real array parallelism (~90s/file PER TASK, but 60 tasks run at
# once under -tc 60). A full scan was already shown to be barely slower than subsampling would be
# (decompression, not parsing, is the bottleneck) -- no reason to subsample now that this isn't
# competing for one core.
#
# Output: results/qc/data/per_sample_global_methylation/<sample>.tsv (sample, hap, n_meth,
# n_unmeth, n_reads) -- one skip-if-exists file per donor, aggregated by the notebook afterward.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
MANIFEST=$PROJDIR/tsv/meta/hprc2_sample_manifest.tsv
MODBED_DIR=$PROJDIR/data/modbed
OUTDIR=$PROJDIR/results/qc/data/per_sample_global_methylation

mkdir -p "$OUTDIR" "$PROJDIR/logs"

# ID=1  # hardcoded for testing; production uses SGE_TASK_ID
ID=$SGE_TASK_ID

ROW=$(awk -F'\t' -v t="$ID" 'NR==t+1 {print}' "$MANIFEST")
if [ -z "$ROW" ]; then
    echo "ERROR: no manifest row for task $ID"
    exit 1
fi
SAMPLE=$(echo "$ROW" | cut -f1)
OUT=$OUTDIR/${SAMPLE}.tsv

echo "$(date): M01 task $ID — $SAMPLE"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

AWK_SCRIPT='{nm=($7==""||$7=="."||$7=="-")?0:gsub(",",",",$7)+1; nu=($8==""||$8=="."||$8=="-")?0:gsub(",",",",$8)+1; tm+=nm; tu+=nu; nr++} END{print tm"\t"tu"\t"nr}'

echo -e "sample\thap\tn_meth\tn_unmeth\tn_reads" > "${OUT}.partial"
for HAP in 1 2; do
    MODBED=$MODBED_DIR/${SAMPLE}_hap${HAP}.modbed.gz
    if [ ! -s "$MODBED" ]; then
        echo "$(date): WARNING no modbed for $SAMPLE hap$HAP, skipping that hap"
        continue
    fi
    time COUNTS=$(zcat "$MODBED" | awk -F'\t' "$AWK_SCRIPT")
    echo -e "${SAMPLE}\t${HAP}\t${COUNTS}" >> "${OUT}.partial"
done
mv "${OUT}.partial" "$OUT"

echo "$(date): $SAMPLE complete"
