#!/bin/bash
set -ue
REFERENCE=${PWD}/result/NC_037304_1_split.fa
qsub -sync yes ./script/run_workflow.job\
     /data/ban-m/a_thaliana/sequel_reads/sequel.fasta\
     ${PWD}/result/pacbio\
     ${REFERENCE}

# /data/ban-m/a_thaliana/genome/mitochondria_enhanced_reference.fa
READ_DIR=/grid/ban-m/arabidopsis_thaliana/sequel/MPI_dataset/
for accession in an1 c24 cvi eri kyo ler sha
do
    qsub -sync yes ./script/run_workflow.job \
         ${READ_DIR}/${accession}.fasta \
         ${PWD}/result/${accession}/ \
         ${REFERENCE}
done


READ_DIR=/grid/ban-m/Arimura/dataset/
for accession in col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    qsub -sync yes ./script/run_workflow.job \
         ${READ_DIR}/${accession}.fasta \
         ${PWD}/result/${accession}/ \
         ${REFERENCE}
done



