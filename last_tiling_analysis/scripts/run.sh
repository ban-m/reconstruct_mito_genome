#!/bin/bash
cargo run --release --bin create_dataset \
      -- ~/data/local_experiment/tiling/canu/mappings.entire.tab \
      ~/data/peaks/guided_asm_sequel_canu.contigs.fasta.mapback.peaks.tsv\
      ~/data/contigs/racon/canu.0.fa\
      ~/data/reads/filtered_reads.fasta
