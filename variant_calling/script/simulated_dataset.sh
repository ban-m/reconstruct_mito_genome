#!/bin/bash


DATA_DIR=../create_dataset_and_circos_plots/result/
ROOT=${PWD}
READ=${DATA_DIR}/easy_simulated_genome/filtered_read.fasta
OUTPATH=${ROOT}/result/easy_simulated_genome/
mkdir -p ${OUTPATH}
REFERENCE=~/work/arabidopsis_mitochondria/create_simulation_data/data/easy/mock_genome_ref.fa
TAB=${DATA_DIR}/easy_simulated_genome/last_db/initial.tab
MAF=${DATA_DIR}/easy_simulated_genome/last_db/initial.maf
qsub -sync yes ./script/workflow_minimal.job \
     ${READ}\
     ${OUTPATH}\
     ${REFERENCE} \
     ${TAB}\
     ${MAF}
