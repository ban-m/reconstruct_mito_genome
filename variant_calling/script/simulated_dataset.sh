#!/bin/bash
DATA_DIR=../create_dataset_and_circos_plots/result/
ROOT=${PWD}
READ=${DATA_DIR}/short_extreme/filtered_read.fasta
OUTPATH=${ROOT}/result/short_extreme/
mkdir -p ${OUTPATH}
REFERENCE=~/work/arabidopsis_mitochondria/create_simulation_data/data/short_extreme/mock_genome_ref.fa
TAB=${DATA_DIR}/short_extreme/last_db/initial.tab
MAF=${DATA_DIR}/short_extreme/last_db/initial.maf
qsub -sync yes ./script/workflow_minimal.job \
     ${READ}\
     ${OUTPATH}\
     ${REFERENCE} \
     ${TAB}\
     ${MAF}
