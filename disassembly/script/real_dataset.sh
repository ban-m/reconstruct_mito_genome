#!/bin/bash
set -ue
REFERENCE=${PWD}/data/NC_037304_1.fa

ROOT=${PWD}
qsub -o ./logfiles/pacbio.log -j y ./script/disassembly.job\
     ${REFERENCE} \
     /data/ban-m/a_thaliana/sequel_reads/sequel.fasta\
     ${PWD}/result/pacbio\
     4 5000

READ_DIR=/grid/ban-m/arabidopsis_thaliana/sequel/MPI_dataset/
for accession in an1 c24 cvi eri kyo ler sha
do
    qsub -o ./logfiles/${accession}.log -j y ./script/disassembly.job \
         ${REFERENCE}\
         ${READ_DIR}/${accession}.fasta \
         ${PWD}/result/${accession}/ \
         4 5000
done

READ_DIR=/grid/ban-m/Arimura/dataset/
for accession in col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    qsub -o ./logfiles/${accession}.log -j y ./script/disassembly.job \
         ${REFERENCE}\
         ${READ_DIR}/${accession}.fasta \
         ${PWD}/result/${accession}/ \
         4 5000
done

