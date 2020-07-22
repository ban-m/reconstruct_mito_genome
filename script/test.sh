#!/bin/bash
set -ue
READ=${PWD}/data/synthetic_data/150_001/reads_150_001.fq
REF=${PWD}/data/synthetic_data/150_001/ref_150_001.fa
mkdir -p ${PWD}/sandbox
BAM=${PWD}/sandbox/temp.bam
CAND=${PWD}/sandbox/cand.vcf
GENO=${PWD}/sandbox/genotype.vcf
PHASED=${PWD}/sandbox/phased.vcf

minimap2 -a -t 23 -x map-pb ${REF} ${READ} | samtools sort -O BAM -m2G -@3 > ${BAM}
samtools index ${BAM}
whatshap find_snv_candidates --multi-allelics -o ${CAND} --pacbio ${REF} ${BAM}
whatshap genotype --ignore-read-groups --reference ${REF} ${CAND} ${BAM} > ${GENO}
whatshap phase --ignore-read-groups --indels --reference ${REF} -o ${PHASED} ${GENO} ${BAM}

