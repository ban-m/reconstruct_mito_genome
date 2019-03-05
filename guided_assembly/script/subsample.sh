#!/bin/bash

### ==== Subsampling BAM file ======

ONT_INPUT=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/mapback
SEQUEL_INPUT=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback

ONT_OUTPUT=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/mapback_subsample
SEQUEL_OUTPUT=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback_subsample

rm -rf ${ONT_OUTPUT} ${SEQUEL_OUTPUT}
mkdir -p ${ONT_OUTPUT} ${SEQUEL_OUTPUT}

function subsample () {
    samtools view -bSh -s 0.3 $1 | samtools sort -m 2G -@2 -O BAM > $2
}

export -f subsample

find ${ONT_INPUT} -name "*.bam" | parallel subsample {} ${ONT_OUTPUT}/{/.}_03.bam
find ${SEQUEL_INPUT} -name "*.bam" | parallel subsample {} ${SEQUEL_OUTPUT}/{/.}_03.bam 

tar -czvf ./result/mapback_subsample03.tar.gz ${ONT_OUTPUT}/* ${SEQUEL_OUTPUT}/*
