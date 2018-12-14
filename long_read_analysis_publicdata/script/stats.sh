#!/bin/bash

function samtoolsstats () {
    samtools stats ${1} > ${2}
}

export -f samtoolsstats
find /grid/ban-m/vigna_species/mapped_reads/ -name *.bam | parallel samtoolsstats {} {.}.stats
