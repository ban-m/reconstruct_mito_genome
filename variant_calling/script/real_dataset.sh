#!/bin/bash

DATA_DIR=../create_dataset_and_circos_plots/result/
ROOT=${PWD}
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    READ=${DATA_DIR}/${accession}/filtered_read.fasta
    OUTPATH=${ROOT}/result/${accession}
    mkdir -p ${OUTPATH}
    REFERENCE=${DATA_DIR}/NC_037304_1_split.fa
    TAB=${DATA_DIR}/${accession}/last_db/initial.tab
    MAF=${DATA_DIR}/${accession}/last_db/initial.maf
    qsub -sync yes ./script/workflow_minimal.job \
         ${READ}\
         ${OUTPATH}\
         ${REFERENCE} \
         ${TAB}\
         ${MAF}
done
