#!/bin/bash
LASTTAB=/grid/ban-m/arabidopsis_thaliana/sequel/tiling/canu/mappings.entire.tab
READ=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fasta
CONTIGS=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/racon/canu.0.fa
cargo run --release --bin main -- ~/data/local_experiment/tiling/canu/mappings.entire.tab \
      ~/data/contigs/racon/canu.0.fa\
      ~/data/reads/filtered_reads.fasta

