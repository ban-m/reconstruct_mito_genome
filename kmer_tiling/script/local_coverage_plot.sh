#!/bin/bash
ASSIGNMENT=./result/sequel_minimap2_assignment.tsv
for k in 10 16 17 18 20 30
do
    Rscript --vanilla --slave ./script/local_coverage_plot.R ./result/${k}_exampler.dat ${ASSIGNMENT} ${k}
done
