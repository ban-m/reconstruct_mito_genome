#!/bin/bash
ROOT=/work/ban-m/arabidopsis_mitochondria
qsub ./scripts/mmmm.job \
     ${ROOT}/create_dataset_and_circos_plots/result/pacbio/filtered_read.fasta \
     ${ROOT}/mito_assembler/result/pacbio/ \
     /data/ban-m/a_thaliana/genome/mitochondria_enhanced_reference.fa

     
