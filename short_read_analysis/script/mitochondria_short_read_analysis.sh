#!/bin/bash
#$ -S /bin/bash
#$ -N SMito
#$ -cwd
#$ -pe smp 12
#$ -e ./log
#$ -o ./out
#$ -V
#$ -m e

ILLUMINA=/data/ban-m/a_thaliana/mitochondria/short_read/ERR2721961.sam
OUTDIR=/data/ban-m/a_thaliana/mitochondria/short_read/
MITO=mitochondria_mapped

function filter () {
    READ=$1
    OUTPUT=$2
    SAMPLE_RATE=$3
    samtools view -H ${READ} > ${OUTPUT}.sam
    samtools view -@ 12 ${READ} | grep "mitochondria"  >> ${OUTPUT}.sam
    samtools sort -@ 12 -m 5G -O BAM ${OUTPUT}.sam > ${OUTPUT}.bam
    samtools view -hb -s ${SAMPLE_RATE} ${OUTPUT}.bam > \
             ${OUTPUT}_subsample${SAMPLE_RATE}.bam
    samtools index ${OUTPUT}_subsample${SAMPLE_RATE}.bam
}

filter $ILLUMINA ${OUTDIR}/illumina_a_thaliana_mito 0.1

