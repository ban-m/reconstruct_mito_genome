#!/bin/bash
echo -e "name\ttotal_base\tn_50\tmean\t#_of contigs\t#_of_n\t#_of_gaps" > ./result/assembly_stats.tsv

S_CANU=/grid/ban-m/arabidopsis_thaliana/sequel/assemble/canu/canu_sequel.contigs.fasta
S_WTDBG=/grid/ban-m/arabidopsis_thaliana/sequel/assemble/wtdbg/wtdbg_sequel.ctg.lay.2nd.fa
S_RA=/grid/ban-m/arabidopsis_thaliana/sequel/assemble/ra/ra_sequel_contigs.fa

N_CANU=/grid/ban-m/arabidopsis_thaliana/nanopore/assemble/canu/canu_nanopore.contigs.fasta
N_WTDBG=/grid/ban-m/arabidopsis_thaliana/nanopore/assemble/wtdbg/wtdbg_nanopore.ctg.lay.2nd.fa
N_RA=/grid/ban-m/arabidopsis_thaliana/nanopore/assemble/ra/ra_ont_contigs.fa

cargo run --release -- ${S_CANU} sequel_canu >> ./result/assembly_stats.tsv
cargo run --release -- ${S_WTDBG} sequel_wtdbg >> ./result/assembly_stats.tsv
cargo run --release -- ${S_RA} sequel_ra >> ./result/assembly_stats.tsv

# cargo run --release -- ${N_CANU} nanopore_canu >> ./result/assembly_stats.tsv
cargo run --release -- ${N_WTDBG} nanopore_wtdbg >> ./result/assembly_stats.tsv
cargo run --release -- ${N_RA} nanopore_ra >> ./result/assembly_stats.tsv

