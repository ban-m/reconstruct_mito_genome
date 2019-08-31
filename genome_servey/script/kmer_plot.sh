#!/bin/bash

for k in $(seq 10 20)
do
    Rscript --vanilla --slave ./script/kmer_plot.R ./result/${k}_mito.hist ${k}_mito
    Rscript --vanilla --slave ./script/kmer_plot.R ./result/${k}_all.hist ${k}_all
done
