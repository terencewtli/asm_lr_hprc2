#!/bin/bash
# M02_postliftover_hap_balance_array.sh -- per-donor hap1-vs-hap2 post-liftover CpG-count/support
# balance check (QC09's cross-check against chain-gap asymmetry). Converted from an in-shell
# sequential loop (single-core, ~2min/donor) to a proper array so it (a) parallelizes and
# (b) survives this interactive session expiring -- see JOURNAL.md 2026-09-16 entry on
# session-persistence risk.
#
# One task per donor in txt/samples/postliftover_full_p03_donors.txt (donors with complete
# 22-chromosome P03 output as of when that list was built -- does not grow automatically as
# P03 progresses further; rebuild the list and resubmit for additional donors later if wanted).
#
# Output: results/qc/data/postliftover_per_donor/<sample>.tsv (one row: sample, hap1_cpgs,
# hap2_cpgs, hap1_support, hap2_support) -- merge afterward with:
#   (echo -e "sample\thap1_cpgs\thap2_cpgs\thap1_support\thap2_support"; \
#    cat results/qc/data/postliftover_per_donor/*.tsv) > results/qc/data/postliftover_hap_balance_all95.tsv

#$ -N M02_postliftover_hap_balance
#$ -cwd
#$ -l h_data=2G,h_rt=0:30:00
#$ -pe shared 1
#$ -t 1-95:1
#$ -tc 60
#$ -o /u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2/logs/M02_postliftover_hap_balance.$JOB_ID.$TASK_ID
#$ -j y

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/asm_lr_hprc2
P03DIR=$PROJDIR/results/all_donors/per_sample_chrom
DONORS=$PROJDIR/txt/samples/postliftover_full_p03_donors.txt
OUTDIR=$PROJDIR/results/qc/data/postliftover_per_donor

sample=$(sed -n "${SGE_TASK_ID}p" "$DONORS")
OUT=$OUTDIR/${sample}.tsv

if [ -s "$OUT" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

echo "$(date): $sample"
zcat ${P03DIR}/${sample}_chr*_cpg_matrix.tsv.gz 2>/dev/null | awk -v s="$sample" '
    NR==1 && $1=="sample" {next}
    $2==1 {c1++; sup1+=($5+$6)}
    $2==2 {c2++; sup2+=($5+$6)}
    END {print s"\t"c1"\t"c2"\t"sup1"\t"sup2}
' > "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"

echo "$(date): $sample done"
