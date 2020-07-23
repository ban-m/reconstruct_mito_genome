#!/bin/bash
# bash ./filter_reads.sh ${READ} ${REFERENCE} > ${FILTERED}
PATH="${PWD}/target/release/:${PATH}"
READ=$1
REFERENCE=$2
minimap2 -t 24 -x map-pb ${REFERENCE} ${READ} |\
    select_mito_reads ${READ} |\
    filter_low_quality |\
    clip_self_chimera 
