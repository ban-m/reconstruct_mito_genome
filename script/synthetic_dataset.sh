#!/bin/bash
ROOT=${PWD}
DATA=${PWD}/data/synthetic_data/
OUTPUT=${PWD}/result/synthetic_data/
mkdir -p ${OUTPUT}
for coverage in 50 100 150
do
    bash ./script/disassembly.sh \
         ${DATA}/${coverage}_001/${coverage}_001_reference.fa \
         ${DATA}/${coverage}_001/${coverage}_001_reads.fa \
         ${OUTPUT}/${coverage}_001/ 2 3000 6 \
         2> ${OUTPUT}/{coverage}_001.log
done
                               

