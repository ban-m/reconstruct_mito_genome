#!/bin/bash
LASTTAB=/grid/ban-m/arabidopsis_thaliana/sequel/tiling/canu/mappings.entire.tab
READ=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fasta
CONTIGS=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/racon/canu.0.fa
# cargo run --release --bin main -- ${LASTTAB} ${CONTIGS} ${READ} 2> ./logfiles/test.log  > ./logfiles/test.out
cargo run --release --bin encode -- ${LASTTAB} ${CONTIGS} ${READ} ./data/contigs.json ./data/reads.json \
      2> ./logfiles/test.log  > ./logfiles/test.out
