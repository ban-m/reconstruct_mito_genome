#!/bin/bash
REF=/data/ban-m/a_thaliana/genome/mitochondria_enhanced_reference.fa
cargo run --release --bin main -- ${REF}
mkdir -p ./data/flip_repeats/
cat ./data/flip_repeats/forward_repeat.fasta |\
    paste - - |\
    grep "original" |\
    tr '\t' '\n' > ./data/flip_repeats/reference.fasta

badread simulate \
        --reference ./data/flip_repeats/forward_repeat.fasta \
        --quantity 500x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --length 15000,5000 > ./data/flip_repeats/forward_repeat.fastq
badread simulate \
        --reference ./data/flip_repeats/reverse_repeat.fasta \
        --quantity 500x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --length 15000,5000 > ./data/flip_repeats/reverse_repeat.fastq

cat ./data/flip_repeats/forward_repeat.fastq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/flip_repeats/read_forward_repeat.fa

cat ./data/flip_repeats/reverse_repeat.fastq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/flip_repeats/read_reverse_repeat.fa

