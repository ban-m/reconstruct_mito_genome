#!/bin/bash

ROOT=/work/ban-m/arabidopsis_mitochondria/create_simulation_data/data
OUTPUT=${PWD}/result/
qsub -o ./logfiles/short_hard.log -j y ./script/disassembly.job \
     ${ROOT}/short_hard/mock_genome_ref.fa \
     ${ROOT}/short_hard/reads.fa \
     ${OUTPUT}/short_hard/ 2

qsub -o ./logfiles/middle_hard.log -j y ./script/disassembly.job \
     ${ROOT}/middle_hard/mock_genome_ref.fa \
     ${ROOT}/middle_hard/reads.fa \
     ${OUTPUT}/middle_hard/ 2

qsub -o ./logfiles/long_hard.log -j y ./script/disassembly.job \
     ${ROOT}/long_hard/mock_genome_ref.fa \
     ${ROOT}/long_hard/reads.fa \
     ${OUTPUT}/long_hard/ 2


