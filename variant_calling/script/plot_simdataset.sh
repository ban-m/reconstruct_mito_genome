#!/bin/bash
# Rscript --vanilla --slave \
#         ./script/plots.R result/easy_simulated_genome/minor_allel_freq.tsv easy_simulated_genome > \
#         ./result/easy_simulated_genome/big_maf_counts.txt

Rscript --vanilla --slave \
        ./script/plots.R \
        result/short_extreme/minor_allel_freq.tsv short_extreme \
        ./result/short_extreme/big_maf_counts.txt
