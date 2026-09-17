#!/bin/bash
#$ -cwd
#$ -o logs/A01c_pmds.$JOB_ID.$TASK_ID
#$ -j y
#$ -N A01c_pmds
#$ -l h_data=2G,h_rt=2:00:00
#$ -pe shared 2
#$ -t 1-40:1

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

source ~/.bashrc

conda activate dnmtools_env

PROJDIR=/u/project/cluo/terencew/igvf/2023_YR2/snmCT/pseudobulk_allc
cd $PROJDIR

# ID=1
ID=${SGE_TASK_ID}

ANNOTS=$PROJDIR/txt/allc_table/rna_cluster_donor.txt
ANNOT=$(head -${ID} $ANNOTS | tail -1)

# METH=$PROJDIR/methcounts/sym_cpg_autosome/rna_cluster_donor/${ANNOT}.methcounts.tsv.gz
METH=$PROJDIR/methcounts/cpg/rna_cluster_donor/${ANNOT}.methcounts.tsv.gz

OUTDIR=$PROJDIR/pmds/dnmtools/rna_cluster_donor/${ANNOT}
mkdir -p $OUTDIR

SUMMARY=$OUTDIR/${ANNOT}.summary
POST=$OUTDIR/${ANNOT}.posteriors
PARAMS=$OUTDIR/${ANNOT}.params

OUT=$OUTDIR/${ANNOT}.bed
SEED=1

time dnmtools pmd -o $OUT -s $SEED -v -i 1000 \
    -S $SUMMARY -r $POST -p $PARAMS $METH

###
PMDS=$PROJDIR/pmds/dnmtools/rna_cluster_donor/${ANNOT}/${ANNOT}.bed
OUT=$PROJDIR/pmds/dnmtools/rna_cluster_donor/${ANNOT}/${ANNOT}.gcov

CHROM_SIZES=/u/project/cluo/terencew/reference/hg38_igvf/GRCh38.chromsizes
time bedtools genomecov -i $PMDS -g $CHROM_SIZES > $OUT

echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "
