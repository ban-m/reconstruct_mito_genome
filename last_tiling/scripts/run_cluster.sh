#!/bin/bash
LASTTAB=/grid/ban-m/arabidopsis_thaliana/sequel/tiling/canu/mappings.entire.tab
READ=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fasta
CONTIGS=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/racon/canu.0.fa
SELF_LASTTAB=/grid/ban-m/arabidopsis_thaliana/sequel/tiling/canu/self.tab
cargo run --release --bin annotate_dotplot -- ${CONTIGS} ${SELF_LASTTAB} > ./data/repeats.json
cargo run --release --bin encode -- ${LASTTAB} ${CONTIGS} ${READ} ./data/contigs.json ./data/reads.json \
      2> ./logfiles/test.log  > ./logfiles/test.out
cargo run --release --bin convert_to_d3_data -- ./data/contigs.json ./data/reads.json > ./data/d3.json

