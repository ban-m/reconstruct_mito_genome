#!/bin/bash
for accession in an1 c24 cvi eri kyo ler sha pacbio col0_1106_exp2
do
    echo ${accession}
    read=${PWD}/result/${accession}/filtered_read/filtered_read.fa
    bbb sequence_statistics --input ${read} --format tsv
done

