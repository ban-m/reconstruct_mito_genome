#!/bin/bash

ROOT=/work/ban-m/arabidopsis_mitochondria/create_simulation_data/data
OUTPUT=${PWD}/result/
qsub -o ./logfiles/short_extreme.log -j y ./script/disassembly.job  \
     ${ROOT}/short_extreme/mock_genome_ref.fa \
     ${ROOT}/short_extreme/reads.fa \
     ${OUTPUT}/short_extreme/ 2 

qsub -o ./logfiles/middle_extreme.log -j y ./script/disassembly.job \
     ${ROOT}/middle_extreme/mock_genome_ref.fa \
     ${ROOT}/middle_extreme/reads.fa \
     ${OUTPUT}/middle_extreme/ 2 

qsub -o ./logfiles/long_extreme.log -j y ./script/disassembly.job \
     ${ROOT}/long_extreme/mock_genome_ref.fa \
     ${ROOT}/long_extreme/reads.fa \
     ${OUTPUT}/long_extreme/ 2 


