#!/bin/bash
set -ue

OUTPUT=${PWD}/data/synthetic_data/
coverage=150
cat ${OUTPUT}/${coverage}_001/reads_${coverage}_001.fq |\
    paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' \
                          > ${OUTPUT}/${coverage}_001/reads_${coverage}_001.fa
