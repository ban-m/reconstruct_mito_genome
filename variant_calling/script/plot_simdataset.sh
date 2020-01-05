#!/bin/bash
# Rscript --vanilla --slave \
#         ./script/plots.R result/easy_simulated_genome/minor_allel_freq.tsv easy_simulated_genome > \
#         ./result/easy_simulated_genome/big_maf_counts.txt

Rscript --vanilla --slave \
        ./script/plots.R result/simulated_genome_long_hard/minor_allel_freq.tsv simulated_genome_long_hard \
        ./result/simulated_genome_long_hard/big_maf_counts.txt
