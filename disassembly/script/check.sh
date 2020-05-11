#!/bin/bash
set -ue
for accession in an1 c24 cvi eri kyo ler sha col0_1106_exp2
do
    input=${PWD}/result/${accession}/multipartite.fasta
    result=$(bbb sequence_statistics --input ${input} --format tsv)
    echo -e "${accession}\n${result}"
done
