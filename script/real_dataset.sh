#!/bin/bash
set -ue
REFERENCE=${PWD}/data/NC_037304_1.fa
READ_DIR=${PWD}/data/filtered_reads/
ROOT=${PWD}
CORES=$1
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2
do
    bash ./script/disassembly.job\
         ${REFERENCE} \
         ${READ_DIR}/${accession}.fasta\
         ${PWD}/result/${accession}\
         3 500 ${CORES} 2> ${PWD}/result/${accession}.log
done
