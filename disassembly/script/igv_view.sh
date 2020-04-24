#$ -S /bin/bash
#$ -N Aln
#$ -cwd
#$ -pe smp 4
#$ -o ./logfiles/alignmnets.log
#$ -j yes
#$ -V
REFERENCE=${PWD}/data/NC_037304_1.fa
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
do
    read=${PWD}/result/${accession}/filtered_read/filtered_read.fa
    minimap2 -t 4 -a -x map-pb ${REFERENCE} ${read} | \
        samtools sort -m 5G -@ 2 -O BAM > ${PWD}/data/${accession}.bam
    samtools index ${PWD}/data/${accession}.bam
done
