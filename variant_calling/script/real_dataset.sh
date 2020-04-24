#$ -S /bin/bash
#$ -N VC
#$ -cwd
#$ -pe smp 1
#$ -o ./logfiles/vc.log
#$ -j yes
#$ -V
#!/bin/bash
set -ue
DATA_DIR=${PWD}/../disassembly/result
GFF_CONVERT=${PWD}/data/gff_convert.tsv
GFF=/grid/ban-m/arabidopsis_thaliana/genome/GCA_000001735.2_TAIR10.1_genomic.gff
cargo build --release
OUTPATH=${PWD}/result/variant_call.tsv
echo "" > ${OUTPATH}
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    READ=${DATA_DIR}/${accession}/filtered_read/filtered_read.fa
    REFERENCE=${PWD}/../disassembly/data/NC_037304_1.fa
    BAM=${PWD}/../disassembly/data/${accession}.bam
    ${PWD}/target/release/variant_calling_longread --bam ${BAM} --reads ${READ} \
          --reference ${REFERENCE} \
          --gff ${GFF} --name ${GFF_CONVERT} -vv |\
        awk -v accession=${accession} 'BEGIN{OFS="\t"}{print $0,accession}' >> ${OUTPATH}
done
