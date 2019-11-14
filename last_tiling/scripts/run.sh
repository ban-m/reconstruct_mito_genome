#!/bin/bash
#$ -S /bin/bash
#$ -N Convert
#$ -cwd
#$ -pe smp 12
#$ -e ./logfiles/convert.log
#$ -o ./logfiles/convert.out
#$ -V
set -ue

# ---- Inputs ---- 
READ=$1
CONTIGS=$2
LASTTAB=$3
SELF_LASTTAB=$4
PREFIX=$5

# ---- Outputs -----

cargo run --release --bin annotate_dotplot -- ${CONTIGS} ${SELF_LASTTAB} \
      > ./data/${PREFIX}_repeats.json
cargo run --release --bin encode \
      -- ${LASTTAB} ${CONTIGS} ${READ} ./data/${PREFIX}_repeats.json \
      ./data/${PREFIX}_contig.json ./data/${PREFIX}_reads.json 
cargo run --release --bin convert_to_d3_data \
      -- ./data/${PREFIX}_contig.json ./data/${PREFIX}_reads.json \
      > ./data/${PREFIX}_d3.json
cat ./scripts/template.html | sed -e "s/SED/${PREFIX}/g" > ./data/${PREFIX}.html
