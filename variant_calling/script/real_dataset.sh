#$ -S /bin/bash
#$ -N VC
#$ -cwd
#$ -pe smp 1
#$ -o ./logfiles/vc.log
#$ -j yes
#$ -V
#!/bin/bash
set -ue
DATA=${PWD}/../disassembly
REFERENCE=${DATA}/data/NC_037304_1.fa
OUTPUT=${PWD}/result/real_dataset.tsv
rm ${OUTPUT}

TAB=${DATA}/result/pacbio/last_db/alignments.tab
READ=${DATA}/result/pacbio/filtered_read/filtered_read.fa
cargo run --release --bin calc_error_rate --\
      ${TAB} ${READ} ${REFERENCE} ler_pacbio >> ${OUTPUT}

for accession in an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    TAB=${DATA}/result/${accession}/last_db/alignments.tab
    READ=${DATA}/result/${accession}/filtered_read/filtered_read.fa
    cargo run --release --bin calc_error_rate --\
          ${TAB} ${READ} ${REFERENCE} ${accession} >> ${OUTPUT}
    
done
