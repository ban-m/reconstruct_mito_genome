#!/bin/bash
DATA=/home/ban-m/work/arabidopsis_mitochondria/genome_servey/data/mitochondria.fa

echo -e 'frequency\tcount\tkmer_size' > ./result/repeat_result.tsv
for mer in 25 50 75 100 125 150 175 200 250 300 350 400 450 500 600 1000 1200
do
    jellyfish count -m ${mer} -s 100M -t 10 -C ${DATA}
    jellyfish histo mer_counts.jf | sed -e "s/$/ ${mer}/g" | sed -e 's/ /\t/g' >> ./result/repeat_result.tsv
done
cat ./result/intersparsed_repeat.maf | awk '{if ( ~ /^a/){print -bash}else{print ,,,,}}'
