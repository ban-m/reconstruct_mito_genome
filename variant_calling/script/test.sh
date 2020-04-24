#!/bin/bash
set -ue
OUTPATH=${PWD}/result/variant_call.tsv
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    echo ${accession}
    cat ${OUTPATH} | rg ${accession} | awk -v accession=${accession} '($5 != "-")'  | wc -l
done
