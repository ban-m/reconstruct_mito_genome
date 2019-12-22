#!/bin/bash
set -ue

LEN=20000
mkdir -p ./data/easy/
cargo run --release --bin create_mock_genomes -- ${LEN} > ./data/easy/mock_genome.fa
badread simulate \
        --reference ./data/easy/mock_genome.fa \
        --quantity 100x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,1000 > ./data/easy/mock_genome_read.fq
cat ./data/easy/mock_genome_read.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/easy/mock_genome_read.fa
cat ./data/easy/mock_genome.fa | paste - - | head -n1 | tr '\t' '\n' > ./data/easy/mock_genome_ref.fa

exit 0

LEN=200000
cargo run --release --bin create_mock_genomes -- ${LEN} > ./data/mock_genome.fa
badread simulate \
        --reference ./data/mock_genome.fa \
        --quantity 50x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ./data/mock_genome_read.fq
cat ./data/mock_genome_read.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/mock_genome_read.fa
cat ./data/mock_genome.fa | paste - - | head -n1 | tr '\t' '\n' > ./data/mock_genome_ref.fa

cargo run --release --bin create_complex_structures -- ${LEN} ./data/complex/
badread simulate \
        --reference ./data/complex/complex1.fa \
        --quantity 100x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ./data/complex/read_complex1.fq
cat ./data/complex/read_complex1.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/complex/read_complex1.fa

badread simulate \
        --reference ./data/complex/complex2.fa \
        --quantity 100x --error_model pacbio \
        --qscore_model pacbio --identity 85,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ./data/complex/read_complex2.fq
cat ./data/complex/read_complex2.fq | paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' > ./data/complex/read_complex2.fa
