#!/bin/bash
#$ -N A02a_methcounts_to_bigwig
#$ -cwd
#$ -l h_data=4G,h_rt=1:00:00
#$ -pe shared 4
#$ -t 1-458:1
#$ -tc 60
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/A02a_methcounts_to_bigwig.$JOB_ID.$TASK_ID
#$ -j y

# Per-(sample,hap) hg38 methylation bigWig from A01a's symmetrized methcounts (NOT P03's
# het-filtered matrix -- see A02a_methcounts_to_bigwig.py's docstring for why). Same 229x2=458
# indexing as A01a_modbed_to_methcounts_array.sh (hold removed 2026-09-17: A01a is complete), so this
# can be submitted any time and will simply wait.
#
# Validated 2026-09-16 on NA19338 (87 real PMDs on chr20 per dnmtools, mean chr20 methylation
# 47.7%, window std 5.9) vs HG01981 (0 PMDs, mean 53.3%, window std 4.8) -- real, sensible,
# donor-differentiated signal, not a flat/garbage track. liftOver mapping rate ~66-70% of
# positions per haplotype (consistent with this session's chain-gap SV-scale-gap measurements,
# ~5-6% of genome per haplotype in gapped regions, plus routine small-indel loss).
#
# Output: data/pmds/<sample>/<sample>_hap{1,2}.hg38.{meth,depth}.bw
# 2026-09-17: rewritten to use the vectorized chain mapping (liftOver version timed out at 1h on
# every task, job 14772524). The 6 old liftOver-based *.hg38.bw files are superseded.

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

OUTDIR=$PROJDIR/data/pmds/$SAMPLE
OUT=$OUTDIR/${SAMPLE}_hap${HAP}.hg38.meth.bw
OUT_DEPTH=$OUTDIR/${SAMPLE}_hap${HAP}.hg38.depth.bw

echo "$(date): A02a task $ID — $SAMPLE hap$HAP"

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

cd "$PROJDIR/scripts/bigwig"
time python3 A02a_methcounts_to_bigwig.py "$SAMPLE" "$HAP" "$OUT" "$OUT_DEPTH"

echo "$(date): $SAMPLE hap$HAP complete"
