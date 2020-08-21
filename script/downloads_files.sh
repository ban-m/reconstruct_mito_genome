#!/bin/bash
ROOT=${PWD}
mkdir -p ${PWD}/data/filtered_reads/
cd ${PWD}/data/filtered_reads/
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2
do
    wget https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/filtered_reads/${accession}.fasta
done
cd ../
wget http://togows.dbcls.jp/entry/nucleotide/NC_037304.1.fasta
cd ../
