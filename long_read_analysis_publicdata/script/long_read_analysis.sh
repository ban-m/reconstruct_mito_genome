#!/bin/bash
#$ -S /bin/bash
#$ -N LRA
#$ -cwd
#$ -pe smp 12
#$ -e ./logfiles/long_read_analysis.log
#$ -o ./logfiles/long_read_analysis.out
#$ -V
#$ -m e

set -ue 
## DATA
SEQUEL=/grid/ban-m/arabidopsis_thaliana/sequel/mapping/sequel_minialign.sam
ONT=/grid/ban-m/arabidopsis_thaliana/nanopore/mapping/nanopore_minialign.sam

function procedure () {
    READ=$1
    OUTPUT=$2
    SAMPLE_RATE=$3
    samtools view -bhS -@12 ${READ} | \
        samtools sort -@12 -m 5G -O BAM > ${OUTPUT}.bam
    samtools index ${OUTPUT}.bam
    samtools view -h -O SAM ${OUTPUT}.bam mitochondria > ${OUTPUT}_mito.sam
    samtools view -hb -s ${SAMPLE_RATE} ${OUTPUT}_mito.sam > \
             ${OUTPUT}_mito_subsample${SAMPLE_RATE}.bam
    samtools index ${OUTPUT}_mito_subsample${SAMPLE_RATE}.bam
    samtools depth -a ${OUTPUT}.bam -r mitochondria > \
             ${OUTPUT}_coverage_mito_MAPQ0.wig
    samtools depth -a ${OUTPUT}.bam -r 1 > ${OUTPUT}_coverage_1_MAPQ0.wig
    samtools depth -a -Q 10 ${OUTPUT}.bam -r mitochondria >\
             ${OUTPUT}_coverage_mito_MAPQ10.wig
    samtools depth -a -Q 10 ${OUTPUT}.bam -r 1 > \
             ${OUTPUT}_coverage_1_MAPQ10.wig
}

mkdir -p /grid/ban-m/arabidopsis_thaliana/nanopore/coverage/
mkdir -p /grid/ban-m/arabidopsis_thaliana/sequel/coverage/

procedure $SEQUEL /grid/ban-m/arabidopsis_thaliana/sequel/coverage/sequel_minialign 0.05
procedure $ONT /grid/ban-m/arabidopsis_thaliana/nanopore/coverage/nanopore_minialign 0.3

