#!/bin/bash
set -ue

LEN=500000
OUTPUT=${PWD}/data/synthetic_data/
mkdir -p ${OUTPUT}/mut01
cargo run --release --bin create_complex_structures -- ${LEN} ${OUTPUT}/mut01 0.001
badread simulate \
        --reference ${OUTPUT}/mut01/complex.fa \
        --quantity 200x --error_model pacbio \
        --seed 10\
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ${OUTPUT}/mut01/read_01.fq
cat ${OUTPUT}/mut01/read_01.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ${OUTPUT}/mut01/read_01.fa

mkdir -p ${OUTPUT}/mut02
cargo run --release --bin create_complex_structures -- ${LEN} ${OUTPUT}/mut02 0.002
badread simulate \
        --reference ${OUTPUT}/mut02/complex.fa \
        --quantity 200x --error_model pacbio \
        --seed 10\
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ${OUTPUT}/mut02/read_02.fq
cat ${OUTPUT}/mut02/read_02.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ${OUTPUT}/mut02/read_02.fa


for coverage in 150 100 50
do
    ## Mut 01
    mkdir -p ${OUTPUT}/${coverage}_01
    cargo run --release --bin create_mock_genomes -- ${LEN} 0.01 2 13223 > \
          ${OUTPUT}/${coverage}_01/contigs_${coverage}_01.fa
    cat ${OUTPUT}/${coverage}_01/contigs_${coverage}_01.fa |\
        head -n2 > ${OUTPUT}/${coverage}_01/ref_${coverage}_01.fa
    badread simulate \
            --reference ${OUTPUT}/${coverage}_01/contigs_${coverage}_01.fa \
            --quantity ${coverage}x --error_model pacbio \
            --seed 10\
            --qscore_model pacbio --identity 85,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,5000 > ${OUTPUT}/${coverage}_01/reads_${coverage}_01.fq
    cat ${OUTPUT}/${coverage}_001/reads_${coverage}_01.fq |\
        paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' \
                              > ${OUTPUT}/${coverage}_01/reads_${coverage}_01.fa
    ## Mut 005
    mkdir -p ${OUTPUT}/${coverage}_005
    cargo run --release --bin create_mock_genomes -- ${LEN} 0.005 2 132121 > \
          ${OUTPUT}/${coverage}_005/contigs_${coverage}_005.fa
    cat ${OUTPUT}/${coverage}_005/contigs_${coverage}_005.fa |\
        head -n2 > ${OUTPUT}/${coverage}_005/ref_${coverage}_005.fa
    badread simulate \
            --reference ${OUTPUT}/${coverage}_005/contigs_${coverage}_005.fa \
            --quantity ${coverage}x --error_model pacbio \
            --seed 10\
            --qscore_model pacbio --identity 85,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,5000 > ${OUTPUT}/${coverage}_005/reads_${coverage}_005.fq
    cat ${OUTPUT}/${coverage}_005/reads_${coverage}_005.fq |\
        paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' \
                              > ${OUTPUT}/${coverage}_005/reads_${coverage}_005.fa
    ## Mut 001
    mkdir -p ${OUTPUT}/${coverage}_001
    cargo run --release --bin create_mock_genomes -- ${LEN} 0.001 2 1321 > \
          ${OUTPUT}/${coverage}_001/contigs_${coverage}_001.fa
    cat ${OUTPUT}/${coverage}_001/contigs_${coverage}_001.fa |\
        head -n2 > ${OUTPUT}/${coverage}_001/ref_${coverage}_001.fa
    badread simulate \
            --reference ${OUTPUT}/${coverage}_001/contigs_${coverage}_001.fa \
            --quantity ${coverage}x --error_model pacbio \
            --seed 10\
            --qscore_model pacbio --identity 85,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,5000 > ${OUTPUT}/${coverage}_001/reads_${coverage}_001.fq
    cat ${OUTPUT}/${coverage}_001/reads_${coverage}_001.fq |\
        paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' \
                              > ${OUTPUT}/${coverage}_001/reads_${coverage}_001.fa
done
