#!/bin/bash
#$ -S /bin/bash
#$ -N Name
#$ -cwd
#$ -pe smp 12
#$ -e ./log
#$ -o ./out
#$ -V
#$ -m e

## DATA
SEQUEL=/data/ban-m/a_thaliana/mitochondria/long_read/a_thaliana_pac.sam
ONT=/data/ban-m/a_thaliana/mitochondria/long_read/a_thaliana_ont.sam
OUTPUT_DIR=/data/ban-m/a_thaliana/mitochondria/long_read
SAMPLE_RATE=0.05_0.3
## Extract only mitochondira mapped reads
function filter () {
    READ=$1
    OUTPUT=$2
    SAMPLE_RATE=$3
    # samtools view -H ${READ} > ${OUTPUT}.sam
    # samtools view -@ 12 ${READ} | grep "mitochondria"  >> ${OUTPUT}.sam
    # samtools sort -@ 12 -m 5G -O BAM ${OUTPUT}.sam > ${OUTPUT}.bam
    samtools view -hb -s ${SAMPLE_RATE} ${OUTPUT}.bam > \
             ${OUTPUT}_subsample${SAMPLE_RATE}.bam
    samtools index ${OUTPUT}_subsample${SAMPLE_RATE}.bam
}

filter $SEQUEL ${OUTPUT_DIR}/sequel_a_thaliana_mito 0.05
filter $ONT ${OUTPUT_DIR}/ont_a_thaliana_mito 0.3

cd ${OUTPUT_DIR}
tar -zvcf long_read_mito${SAMPLE_RATE}.tar.gz *_subsample${SAMPLE_RATE}*.ba[im]

SEQUEL_BAM=/data/ban-m/a_thaliana/read_alignment/a_thaliana_pac.bam
ONT_BAM=/data/ban-m/a_thaliana/read_alignment/a_thaliana_ont.bam

samtools depth -a $SEQUEL_BAM > ${OUTPUT_DIR}/sequel_coverage.wig
samtools depth -a ${ONT_BAM} > ${OUTPUT_DIR}/ont_coverage.wig
