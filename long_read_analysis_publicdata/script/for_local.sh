#!/bin/bash
OUTPUT_DIR=/data/ban-m/a_thaliana/mitochondria/long_read
SEQUEL=sequel_a_thaliana_mito
ONT=ont_a_thaliana_mito
SAMPLE_RATE=0.01
function subsample () {
    OUTPUT=$1
    samtools view -hb -s ${SAMPLE_RATE} ${OUTPUT_DIR}/${OUTPUT}.bam > \
             ${OUTPUT_DIR}/${OUTPUT}.subsample_${SAMPLE_RATE}.bam
    samtools index ${OUTPUT_DIR}/${OUTPUT}.subsample_${SAMPLE_RATE}.bam
}

subsample $SEQUEL
subsample $ONT

cd ${OUTPUT_DIR}
tar -zvcf long_read_mito${SAMPLE_RATE}.tar.gz subsample_${SAMPLE_RATE}*.ba[im]
