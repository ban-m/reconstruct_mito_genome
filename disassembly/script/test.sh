#!/bin/bash
for read in $( find . -maxdepth 3 -size -1M -name "*.fasta" )
do
    echo ${read}
done
