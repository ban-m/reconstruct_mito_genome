#!/bin/bash
set -ue
mkdir -p png pdf data

ONT=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/mapback/
SEQUEL=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback/


function coverage_calc () {
    echo $1
    # samtools depth -d 1000000 $1 > $2
    Rscript --vanilla --slave ./script/coverage.R $2
}

export -f coverage_calc


find ${ONT} -name "*.bam" | parallel coverage_calc {} ./data/{/.}.wig
find ${SEQUEL} -name "*.bam" | parallel coverage_calc {} ./data/{/.}.wig

cat ./data/*.csv  > ./result/coverage_summary.csv
