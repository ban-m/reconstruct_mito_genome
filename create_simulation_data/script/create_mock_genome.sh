#!/bin/bash
set -ue

function create_easy(){
    OUTPATH=$1
    LEN=$2
    mkdir -p ${OUTPATH}
    cargo run --release --bin create_mock_genomes ${LEN} > ${OUTPATH}/mock_genome.fa
    badread simulate \
            --reference ${OUTPATH}/mock_genome.fa \
            --quantity 100x --error_model pacbio \
            --qscore_model pacbio --identity 85,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,1000 > ${OUTPATH}/reads.fq
    cat ${OUTPATH}/reads.fq | paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' > ${OUTPATH}/reads.fa
    cat ${OUTPATH}/mock_genome.fa | paste - - | head -n1 | \
        tr '\t' '\n' > ${OUTPATH}/mock_genome_ref.fa
}

# create_easy ./data/short_easy 20000
# create_easy ./data/middle_easy 200000
# create_easy ./data/long_easy 500000

function create_hard() {
    OUTPATH=$1
    LEN=$2
    mkdir -p ${OUTPATH}
    cargo run --release --bin create_mock_genomes_hard -- ${LEN} > ${OUTPATH}/mock_genome.fa
    badread simulate \
            --reference ${OUTPATH}/mock_genome.fa \
            --quantity 150x --error_model pacbio \
            --qscore_model pacbio --identity 85,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,1000 > ${OUTPATH}/reads.fq
    cat ${OUTPATH}/reads.fq | paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' > ${OUTPATH}/reads.fa
    cat ${OUTPATH}/mock_genome.fa | paste - - | head -n1 | \
        tr '\t' '\n' >  ${OUTPATH}/mock_genome_ref.fa
}

# create_hard ./data/short_hard 20000
# create_hard ./data/middle_hard 200000
# create_hard ./data/long_hard 500000


function create_extreme() {
    OUTPATH=$1
    LEN=$2
    mkdir -p ${OUTPATH}
    cargo run --release --bin create_mock_genomes_extreme -- ${LEN} > ${OUTPATH}/mock_genome.fa
    badread simulate \
            --reference ${OUTPATH}/mock_genome.fa \
            --quantity 200x --error_model pacbio \
            --qscore_model pacbio --identity 85,95,3 \
            --junk_reads 0 --random_reads 0 --chimeras 0 \
            --length 15000,1000 > ${OUTPATH}/reads.fq
    cat ${OUTPATH}/reads.fq | paste - - - - | cut -f 1,2 |\
        sed -e 's/@/>/g' | tr '\t' '\n' > ${OUTPATH}/reads.fa
    cat ${OUTPATH}/mock_genome.fa | paste - - | head -n1 | tr '\t' '\n' > ${OUTPATH}/mock_genome_ref.fa
}

# create_extreme ./data/short_extreme 20000
# create_extreme ./data/middle_extreme 200000
# create_extreme ./data/long_extreme 500000

LEN=400000
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

