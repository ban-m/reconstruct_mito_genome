#!/bin/bash
OUTPUT=${PWD}/result/complex1
READ=${OUTPUT}/filtered_read/filtered_read.fa
REFERENCE=/work/ban-m/arabidopsis_mitochondria/create_simulation_data/data/complex/read_complex1.fa
mmmm create_viewer --assignments ${OUTPUT}/readlist.tsv \
     --contig_aln ${OUTPUT}/allcontigs.aln.tab \
     --contigs ${OUTPUT}/multipartite.fasta \
     --output_dir ${OUTPUT}/viewer/ \
     --read_aln ${OUTPUT}/allreads.aln.tab \
     --reads ${READ} \
     --reference ${REFERENCE}

