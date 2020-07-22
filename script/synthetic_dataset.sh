#!/bin/bash
ROOT=${PWD}
DATA=${PWD}/data/synthetic_data/
OUTPUT=${PWD}/result/synthetic_data/
bash ${PWD}/script/disassembly.job \
     ${DATA}/mut01/reference.fa \
     ${DATA}/mut01/read_01.fa \
     ${OUTPUT}/mut01/ 3 4000 24

bash ${PWD}/script/disassembly.job \
     ${DATA}/mut02/reference.fa \
     ${DATA}/mut02/read_02.fa \
     ${OUTPUT}/mut02/ 3 4000 24

