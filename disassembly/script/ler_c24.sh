#!/bin/bash
set -ue
REFERENCE=/grid/ban-m/arabidopsis_thaliana/genome/JF729202.fasta

ROOT=${PWD}
qsub -o ./logfiles/pacbio_ler.log -j y ./script/disassembly.job\
     ${REFERENCE} \
     /data/ban-m/a_thaliana/sequel_reads/sequel.fasta\
     ${PWD}/result/pacbio_ler\
     4 3600

READ_DIR=/grid/ban-m/arabidopsis_thaliana/sequel/MPI_dataset/
for accession in ler
do
    qsub -o ./logfiles/${accession}_ler_2.log -j y ./script/disassembly.job \
         ${REFERENCE}\
         ${READ_DIR}/${accession}.fasta \
         ${PWD}/result/${accession}_ler_2/ \
         4 3600
done

REFERENCE=/grid/ban-m/arabidopsis_thaliana/genome/JF729200.fasta
for accession in c24
do
    qsub -o ./logfiles/${accession}_c24.log -j y ./script/disassembly.job \
         ${REFERENCE}\
         ${READ_DIR}/${accession}.fasta \
         ${PWD}/result/${accession}_c24/ \
         4 5000
done
