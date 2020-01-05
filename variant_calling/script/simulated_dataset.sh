#!/bin/bash


DATA_DIR=../create_dataset_and_circos_plots/result/
ROOT=${PWD}
READ=${DATA_DIR}/simulated_genome_long_hard/filtered_read.fasta
OUTPATH=${ROOT}/result/simulated_genome_long_hard/
mkdir -p ${OUTPATH}
REFERENCE=~/work/arabidopsis_mitochondria/create_simulation_data/data/long_hard/mock_genome_ref.fa
TAB=${DATA_DIR}/simulated_genome_long_hard/last_db/initial.tab
MAF=${DATA_DIR}/simulated_genome_long_hard/last_db/initial.maf
qsub -sync yes ./script/workflow_minimal.job \
     ${READ}\
     ${OUTPATH}\
     ${REFERENCE} \
     ${TAB}\
     ${MAF}
