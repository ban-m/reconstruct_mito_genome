#!/bin/bash
set -ue
cd /grid/ban-m/arabidopsis_thaliana/sequel/MPI_dataset/
for accession in an1 c24 cvi eri kyo ler sha
do
    cat /grid/ban-m/arabidopsis_thaliana/dataset/sequel/${accession}/* | \
        paste - - - - | \
        cut -f 1,2 | \
        sed -e  's/@/>/g' |\
        tr '\t' '\n' > ${accession}.fasta
done
