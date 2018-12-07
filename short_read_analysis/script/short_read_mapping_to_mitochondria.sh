#!/bin/bash
#$ -S /bin/bash
#$ -N ShortRead
#$ -cwd
#$ -pe smp 12
#$ -e ./logfiles/short_read_alignmnet.log
#$ -o ./logfiles/short_read_alignment.out
#$ -V

READ1=/data/ban-m/a_thaliana/illumina/ERR2721961_1.fastq
READ2=/data/ban-m/a_thaliana/illumina/ERR2721961_2.fastq
GENOME=/data/ban-m/a_thaliana/genome/TAIR10_chr_all.fas
OUTPUT_DIR=/data/ban-m/a_thaliana/mitochondria/short_read/
OUTPUT=ERR2721961.sam

mkdir -p ${OUTPUT_DIR}
bwa index ${GENOME}
bwa mem -t 12 ${GENOME} ${READ1} ${READ2} > ${OUTPUT_DIR}/${OUTPUT}

