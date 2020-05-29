#!/bin/bash
REFERENCE=${PWD}/data/NC_037304_1.fa
READ_DIR=/grid/ban-m/arabidopsis_thaliana/sequel/MPI_dataset/
accession=kyo
qsub -o ./logfiles/${accession}.log -j y ./script/disassembly.job \
     ${REFERENCE}\
     ${READ_DIR}/${accession}.fasta \
     ${PWD}/result/${accession}/ \
     3 5000
