#!/bin/bash

ROOT=/work/ban-m/arabidopsis_mitochondria/create_simulation_data/data
OUTPUT=${PWD}/result/
qsub -o ./logfiles/short_easy.log -j y ./scripts/disassembly.job  \
     ${ROOT}/short_easy/mock_genome_ref.fa \
     ${ROOT}/short_easy/reads.fa \
     ${OUTPUT}/short_easy/ 2 

qsub -o ./logfiles/middle_easy.log -j y ./scripts/disassembly.job \
     ${ROOT}/middle_easy/mock_genome_ref.fa \
     ${ROOT}/middle_easy/reads.fa \
     ${OUTPUT}/middle_easy/ 2 

qsub -o ./logfiles/long_easy.log -j y ./scripts/disassembly.job\
     ${ROOT}/long_easy/mock_genome_ref.fa \
     ${ROOT}/long_easy/reads.fa \
     ${OUTPUT}/long_easy/ 2


