#!/bin/bash
ROOT=${PWD}
DATA=${PWD}/data/synthetic_data/
OUTPUT=${PWD}/result/synthetic_data/
for coverage in 150 100 50
do
    bash ./script/disassembly.sh \
         ${DATA}/${coverage}_001/ref_${coverage}_001.fa \
         ${DATA}/${coverage}_001/reads_${coverage}_001.fa \
         ${OUTPUT}/${coverage}_001/ 3 4000 24
done

