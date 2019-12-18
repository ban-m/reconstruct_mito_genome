#!/bin/bash
#$ -S /bin/bash
#$ -N decompose
#$ -cwd
#$ -pe smp 1
#$ -e ./log
#$ -o ./out
#$ -V
#$ -m e
READ=../create_dataset_and_circos_plots/result/forward_repeat/filtered_read.fasta
ALIGN=../create_dataset_and_circos_plots/result/forward_repeat/last_db/initial.tab
CONTIG=../create_simulation_data/data/flip_repeats/reference.fasta
REPEAT=../create_dataset_and_circos_plots/result/forward_repeat/circos/repeats.json
# cargo run --release --bin test -- ${READ} ${ALIGN} ${CONTIG}
cargo run --release --bin enumerate_cr -- ${READ} ${ALIGN} ${CONTIG} ${REPEAT} \
      > ./logfiles/cr.json
