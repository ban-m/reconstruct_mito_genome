#!/bin/bash
#$ -S /bin/bash
#$ -N Bench
#$ -cwd
#$ -pe smp 24
#$ -V
#$ -e ./logfiles/mitoasm.log
#$ -o ./logfiles/mitoasm.out
#$ -m e
# -M banmasutani@gmail.com ## To send a mail when ended
# -t 1:n ## For array job
set -ue

### Variables

SEQUEL=/work/ban-m/arabidopsis_mitochondria/create_simulation_data/data/forward_repeat.fastq
OUTPUT=${PWD}/result

#### ==== Canu =====
rm -rf ${OUTPUT}/canu/
mkdir -p ${OUTPUT}/canu/
qsub -o ./logfiles/canu_sequel.out -e ./logfiles/canu_sequel.log \
     ./scripts/canu.job \
     ${OUTPUT}/canu \
     canu \
     -pacbio-raw \
     ${SEQUEL}

#### ==== Wtdbg2 =====
rm -rf ${OUTPUT}/wtdbg2
mkdir -p ${OUTPUT}/wtdbg2
qsub -o ./logfiles/wtdbg2_sequel.out -e ./logfiles/wtdbg2_sequel.log \
     ./scripts/wtdbg2.job \
     ${OUTPUT}/wtdbg2 \
     sequel \
     ${SEQUEL} \
     wtdbg2 \
     map-pb 

#### ===== Flye =======
rm -rf ${OUTPUT}/flye
mkdir -p ${OUTPUT}/flye
qsub -o ./logfiles/flye_sequel.out -e ./logfiles/flye_sequel.log \
     ./scripts/flye.job \
     pacbio-raw \
     ${SEQUEL} \
     ${OUTPUT}/flye



