#!/bin/bash
badread simulate \
        --reference ./data/mito_mutant.fasta --quantity 400x \
        --error_model pacbio \
        --qscore_model pacbio \
        --identity 85,95,3 \
        --length 15000,5000 \
        --junk_reads 0 \
        --random_reads 0 \
        --chimeras 0 \
        > ./data/mito_mutant_simulate.fastq

