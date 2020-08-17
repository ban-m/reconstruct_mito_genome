#!/bin/bash
set -ue
ROOT=${PWD}
DATA=${PWD}/data/synthetic_data/
OUTPUT=${PWD}/result/synthetic_data/
coverage=50
for err in 001 
do
    mmmm decompose \
         --output ${OUTPUT}/test/ \
         --reads ${READ} --contigs ${REFERENCE} \
         --cluster_num 1 --threads 12 \
         --limit ${LIMIT}\
         -vv
done
