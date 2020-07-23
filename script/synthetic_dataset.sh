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

for coverage in 150 100 50
do
    bash ./script/disassembly.sh \
         ${DATA}/${coverage}_01/ref_${coverage}_01.fa \
         ${DATA}/${coverage}_01/reads_${coverage}_01.fa \
         ${OUTPUT}/${coverage}_01/ 3 4000 24
    bash ./script/disassembly.sh \
         ${DATA}/${coverage}_005/ref_${coverage}_005.fa \
         ${DATA}/${coverage}_005/reads_${coverage}_005.fa \
         ${OUTPUT}/${coverage}_005/ 3 4000 24
    bash ./script/disassembly.sh \
         ${DATA}/${coverage}_001/ref_${coverage}_001.fa \
         ${DATA}/${coverage}_001/reads_${coverage}_001.fa \
         ${OUTPUT}/${coverage}_001/ 3 4000 24
done

