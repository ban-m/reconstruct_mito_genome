#!/bin/bash
set -ue
ROOT=${PWD}
DATA=${PWD}/data/synthetic_data/
OUTPUT=${PWD}/result/synthetic_data/
coverage=100
qsub -o ./logfiles/${coverage}_001.log -j y -S /bin/bash -cwd -pe smp 24 -V \
     ./script/disassembly.sh \
     ${DATA}/${coverage}_001/${coverage}_001_reference.fa \
     ${DATA}/${coverage}_001/${coverage}_001_reads.fa \
     ${OUTPUT}/${coverage}_001/ 2 4000 24
# REFERENCE=${PWD}/data/NC_037304_1.fa
# READ_DIR=${PWD}/data/filtered_reads/
# ROOT=${PWD}
# mkdir -p logfiles
# accession=pacbio
# OUTPUT=${PWD}/result/pacbio
# READ=${PWD}/data/filtered_reads/pacbio.fasta
# MIN_CLUSTER=3
# CORES=24
# LIMIT=4000
# ${PWD}/target/release/mmmm decompose \
#      --alignments ${OUTPUT}/last_db/alignments.tab --output ${PWD}/result/pacbio_2 \
#      --reads ${READ} --contigs ${REFERENCE} \
#      --self_alignments ${OUTPUT}/last_db/self.tab \
#      --cluster_num ${MIN_CLUSTER} --threads ${CORES} \
#      --limit ${LIMIT} --resume ${OUTPUT}/encoded_reads.json -vv
