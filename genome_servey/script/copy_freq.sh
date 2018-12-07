#!/bin/bash
DATA=/home/ban-m/work/arabidopsis_mitochondria/genome_servey/data/mitochondria.fa

for mer in 50 100 150 200 500
do
    jellyfish count -m ${mer} -s 100M -t 10 -C -o ./result/${mer}_count.jf ${DATA}
    jellyfish histo ./result/${mer}_count.jf > ./result/${mer}_histogram.dat
done
