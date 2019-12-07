#!/bin/bash
#$ -S /bin/bash
#$ -N Convert
#$ -cwd
#$ -pe smp 12
#$ -e ./logfiles/convert.log
#$ -o ./logfiles/convert.out
#$ -V
set -ue
ROOT=/grid/ban-m/arabidopsis_thaliana/sequel/assemble/pacbio_w_reference/
READ=${ROOT}/filtered_read.fasta
REFERENCE=${ROOT}/mito.fa
mkdir -p ./data/reference/
SELF_LASTTAB=${ROOT}/last_db/self.tab
PREFIX=./data/mine_to_reference
LASTTAB=${ROOT}/last_db/initial.tab
cargo run --release --bin main \
      -- ${LASTTAB} ${REFERENCE} ${READ}
