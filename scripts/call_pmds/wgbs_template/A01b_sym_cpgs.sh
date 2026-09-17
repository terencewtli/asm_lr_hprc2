#!/bin/bash
#$ -cwd
#$ -o logs/A01b_sym_cpg.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01b_sym_cpg
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-36:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate allcools

PROJDIR=/u/project/cluo/terencew/igvf/2023_YR2/snmCT/pseudobulk_allc
cd $PROJDIR

# ID=1
ID=${SGE_TASK_ID}

ANNOTS=$PROJDIR/txt/allc_table/rna_cluster_donor.txt
ANNOT=$(head -${ID} $ANNOTS | tail -1)

METH=$PROJDIR/methcounts/cpg/rna_cluster_donor/${ANNOT}.methcounts.tsv.gz
SYM=$PROJDIR/methcounts/sym_cpg/rna_cluster_donor/${ANNOT}.methcounts.tsv.gz

DNMTOOLS=/u/home/t/terencew/project-cluo/programs/dnmtools-1.0.0/dnmtools

time $DNMTOOLS sym -v -o $SYM $METH

AUTO_SYM=$PROJDIR/methcounts/sym_cpg_autosome/rna_cluster_donor/${ANNOT}.methcounts.tsv.gz

time zcat $SYM | grep -v -e 'chrX' -e 'chrY' -e 'chrM' | gzip - > $AUTO_SYM

###
OUT=$PROJDIR/methcounts/levels/rna_cluster_donor/${ANNOT}.levels
time dnmtools levels -relaxed -v -o $OUT $SYM

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
