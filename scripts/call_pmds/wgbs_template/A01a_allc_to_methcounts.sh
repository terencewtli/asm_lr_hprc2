#!/bin/bash
#$ -cwd
#$ -o logs/A01a_allc_to_methcounts.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01a_allc_to_methcounts
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-40:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate allcools

PROJDIR=/u/project/cluo/terencew/igvf/2023_YR2/snmCT/pseudobulk_allc/
cd $PROJDIR

# ID=3
ID=${SGE_TASK_ID}

ANNOTS=$PROJDIR/txt/allc_table/rna_cluster_donor.txt
ANNOT=$(head -${ID} $ANNOTS | tail -1)

ALLC=$PROJDIR/allc/rna_cluster_donor/${ANNOT}.CGN-Merge.allc.tsv.gz
OUT=$PROJDIR/methcounts/cpg/rna_cluster_donor/${ANNOT}.methcounts.tsv.gz

### extract CGN and rename to CpG

time zcat $ALLC | cut -f1-6 - \
    | awk -v OFS='\t' \
    '{ $4 = "CpG" } {$5 /= $6 } { print $0 }' \
    | gzip - > $OUT

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
