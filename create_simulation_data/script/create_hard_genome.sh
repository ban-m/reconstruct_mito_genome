#!/bin/bash
set -ue

function create_mock() {
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
    cat ${OUTPATH}/mock_genome.fa | paste - - | head -n1 | tr '\t' '\n' > ${OUTPATH}/mock_genome_ref.fa
}

create_mock ./data/short_hard 20000
create_mock ./data/middle_hard 200000
create_mock ./data/long_hard 300000
