#!/bin/bash
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    Rscript --vanilla --slave \
            ./script/plots.R result/${accession}/minor_allel_freq.tsv ${accession} > \
            ./result/${accession}/big_maf_counts.txt
done
                 
