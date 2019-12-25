#!/bin/bash
Rscript --vanilla --slave \
        ./script/plots.R result/easy_simulated_genome/minor_allel_freq.tsv easy_simulated_genome > \
        ./result/easy_simulated_genome/big_maf_counts.txt
                 
