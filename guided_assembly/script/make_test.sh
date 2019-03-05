#!/bin/bash
samtools view -s 0.0001 /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback/guided_asm_sequel_canu.contigs.fasta.mapback.bam > ./data/test_canu_sequel.sam
cat ./data/test_canu_sequel.sam | cut -d$'\t' -f -1 > ./name.dat
cat /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fq | \
    paste - - - - | grep -f ./name.dat | tr '\t' '\n'  > \
                                            ./data/test_canu_sequel.fastq

                                           
