#!/bin/bash
reference=${PWD}/data/NC_037304_1.fa
read=${PWD}/data/filtered_reads/pacbio.fasta
file=test
BAM=${PWD}/sandbox/${file}.bam
CAND=${PWD}/sandbox/${file}_cand.vcf
GENOTYPE=${PWD}/sandbox/${file}_genot.vcf
PHASED=${PWD}/sandbox/${file}_phased.vcf.gz
PHASED_BAM=${PWD}/sandbox/${file}.phased.bam
RESULT=${PWD}/sandbox/${file}.tsv
minimap2 -R '@RG\tID:1\tSM:Sample1' -a -t 23 -x map-pb ${reference} ${read} |\
    samtools sort -O BAM -m2G -@3 > ${BAM}
samtools index ${BAM}
whatshap find_snv_candidates --multi-allelics -o ${CAND} --pacbio \
         ${reference}  ${BAM}
whatshap genotype --ignore-read-groups --reference ${reference} \
         ${CAND} ${BAM} > ${GENOTYPE}
whatshap phase --algorithm hapchat --reference ${reference} --ignore-read-groups \
         -o ${PHASED} ${GENOTYPE} ${BAM}
tabix -p vcf ${PHASED}
whatshap haplotag -o ${PHASED_BAM} \
         --ignore-read-groups \
         --output-haplotag-list ${RESULT} \
         --reference ${reference} ${PHASED} ${BAM} \
         
