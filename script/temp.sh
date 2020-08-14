#!/bin/bash
set -ue
PATH=${PWD}/target/release:${PATH}
OUTPUT=${PWD}/result/pacbio_2/
READ=${PWD}/data/filtered_reads/pacbio.fasta
REFERENCE=${PWD}/data/NC_037304_1.fa
MIN_CLUSTER=3
CORES=24
LIMIT=4000
mmmm decompose \
     --output ${OUTPUT} \
     --reads ${READ} --contigs ${REFERENCE} \
     --cluster_num ${MIN_CLUSTER} --threads ${CORES} \
     --limit ${LIMIT} --resume ${PWD}/result/pacbio/encoded_reads.json \
     -vv 2> ${PWD}/logfiles/pacbio_2.log
