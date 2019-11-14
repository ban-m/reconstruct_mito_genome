#!/bin/bash

LEN=300000
cargo run --release --bin create_mock_genomes -- ${LEN} > ./data/mock_genome.fa
badread simulate \
        --reference ./data/mock_genome.fa \
        --quantity 100x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ./data/mock_genome_read.fq
cat ./data/mock_genome_read.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/mock_genome_read.fa
cat ./data/mock_genome.fa | paste - - | head -n1 | tr '\t' '\n' > ./data/mock_genome_ref.fa


