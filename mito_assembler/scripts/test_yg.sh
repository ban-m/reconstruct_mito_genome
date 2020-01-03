#!/bin/bash
set -ue
ROOT=${PWD}
READS=${PWD}/../create_simulation_data/data/mock_genome_read.fa
OUTPATH=${PWD}/result/mock_genome
REFERENCE=${PWD}/../create_simulation_data/data/mock_genome_ref.fa
# qsub ./scripts/yellow_gate.job ${READS} ${OUTPATH} ${REFERENCE}
for reads in ${OUTPATH}/yg/*.fasta
do
    name=${reads%.fasta}
    name=${name##*/}
    # qsub -sync yes \
    #-o ./logfiles/canu_sequel_${name}.out \
    #      -e ./logfiles/canu_sequel_${name}.log \
    #      ./scripts/canu.job \
    #      ${OUTPATH}/canu_${name} \
    #      canu_${name} \
    #      -pacbio-raw \
    #      ${reads}
    ## Correction
    contig=${OUTPATH}/canu_${name}/canu_${name}.contigs.fasta
    mkdir -p ${OUTPATH}/${name}
    qsub -o ./logfiles/racon_${name}.log\
         -e ./logfiles/racon_${name}.log\
         ./scripts/racon.job\
         ${reads} ${contig} ${OUTPATH}/${name}
done

