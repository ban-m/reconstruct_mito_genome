#!/bin/bash
set -ue
LEN=500000
OUTPUT=${PWD}/data/test/
for coverage in 150
do
    mkdir -p ${OUTPUT}/${coverage}_001
    cargo run --release --bin create_mock_genomes -- ${LEN} 0.001 2 13223 \
          ${OUTPUT}/${coverage}_001/${coverage}_001
    badread simulate \
            --reference ${OUTPUT}/${coverage}_001/${coverage}_001_contigs.fa \
            --quantity ${coverage}x --error_model pacbio \
            --seed 10\
            --qscore_model pacbio --identity 90,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,5000 > ${OUTPUT}/${coverage}_001/${coverage}_001_reads.fq
    cat ${OUTPUT}/${coverage}_001/${coverage}_001_reads.fq |\
        paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' \
                              > ${OUTPUT}/${coverage}_001/${coverage}_001_reads.fa
done

DATA=${OUTPUT}
OUTPUT=${PWD}/result/test/
coverage=150
qsub -o ./logfiles/test_${coverage}_001.log -j y -S /bin/bash -cwd -pe smp 23 -V \
     ./script/disassembly.sh \
     ${DATA}/${coverage}_001/${coverage}_001_reference.fa \
     ${DATA}/${coverage}_001/${coverage}_001_reads.fa \
     ${OUTPUT}/${coverage}_001/ 2 3000 23

