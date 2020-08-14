#!/bin/bash
set -ue
ROOT=${PWD}
DATA=${PWD}/data/synthetic_data/
OUTPUT=${PWD}/result/synthetic_data/
coverage=50
for err in 001 
do
    bash \
        ./script/disassembly.sh \
        ${DATA}/${coverage}_${err}/${coverage}_${err}_reference.fa \
        ${DATA}/${coverage}_${err}/${coverage}_${err}_reads.fa \
        ${OUTPUT}/test_${err}/ 2 2000 24 1> ./logfiles/50_001_2.log 2>&1 
done
