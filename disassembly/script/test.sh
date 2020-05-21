#!/bin/bash
#$ -S /bin/bash
#$ -N M4TEST2931
#$ -cwd
#$ -pe smp 24
#$ -V
#$ -m e
#$ -j y
#$ -o ./logfiles/temp_test_2931.log
OUTPUT=${PWD}/result/complex1
TEMPDIR=${PWD}/result/temp2931_complex1
mkdir -p ${TEMPDIR}
mmmm decompose --alignments ${OUTPUT}/last_db/alignments.tab --output ${TEMPDIR} \
     --reads ${OUTPUT}/filtered_read/filtered_read.fa\
     --contigs ${PWD}/../create_simulation_data/data/complex/reference.fa \
     --self_alignments ${OUTPUT}/last_db/self.tab \
     --cluster_num 3 --threads 24 \
     --limit 3000\
     -vv
