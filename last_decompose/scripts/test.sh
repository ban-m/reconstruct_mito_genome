#!/bin/bash
#$ -S /bin/bash
#$ -N decompose
#$ -cwd
#$ -pe smp 1
#$ -e ./log
#$ -o ./out
#$ -V
#$ -m e
ROOT=/grid/ban-m/arabidopsis_thaliana/sequel/assemble/mine
READ=${ROOT}/filtered_read.fasta
ALIGN=${ROOT}/last_db/collupsed.tab
CONTIG=${ROOT}/contigs.fasta
cargo run --release --bin test -- ${READ} ${ALIGN} ${CONTIG}
