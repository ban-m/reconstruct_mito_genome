#!/bin/bash
#$ -S /bin/bash
#$ -N kmer
#$ -cwd
#$ -pe smp 24
#$ -e ./logfiles/kmer-frequency.log
#$ -o ./logfiles/kmer-frequency.out
#$ -V
#$ -m e

# DATA=/home/ban-m/work/arabidopsis_mitochondria/genome_servey/data/mitochondria.fa
# echo -e 'frequency\tcount\tkmer_size' > ./result/repeat_result.tsv
# for mer in 25 50 75 100 125 150 175 200 250 300 350 400 450 500 600 1000 1200
# do
#     jellyfish count -m ${mer} -s 4G --bf-size 50G -t 10 -o ${mer}.jf ${DATA}
#     jellyfish histo ${mer}.jf | sed -e "s/$/ ${mer}/g" | sed -e 's/ /\t/g' >> ./result/repeat_result.tsv
# done
# cat ./result/intersparsed_repeat.maf | awk '{if ( ~ /^a/){print -bash}else{print ,,,,}}'

cd ./result/
GENOME=/grid/ban-m/arabidopsis_thaliana/genome/GCF_000001735.4_TAIR10.1_genomic.fna
MITOCHONDRIA=/grid/ban-m/arabidopsis_thaliana/genome/mitochondria.fa
for k in $(seq 10 20)
do
    jellyfish count -s 4G --bf-size 10G -m ${k} -t 24 -o ${k}_all.jf ${GENOME}
    jellyfish histo -o ${k}_all.hist ${k}_all.jf
    jellyfish count -s 4G --bf-size 10G -m ${k} -t 24 -o ${k}_mito.jf ${MITOCHONDRIA}
    jellyfish histo -o ${k}_mito.hist ${k}_mito.jf
done

