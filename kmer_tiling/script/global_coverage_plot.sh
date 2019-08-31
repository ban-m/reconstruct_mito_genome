#!/bin/bash
ASSIGNMENT=./result/sequel_minimap2_assignment.tsv
for k in 10 16 17 18 20 30
do
    Rscript --vanilla --slave ./script/global_coverage_plot.R ./result/tiling_${k}.tsv ${ASSIGNMENT}
done

