#!/bin/bash


## ===== Files =====

GENBANK=/data/ban-m/mitochondrial_databese/mitochondrion.2.genomic.gbff
ARABIDOPSIS=/data/ban-m/a_thaliana/genome/mitochondria_2line.fa
ARABIDOPSIS_GFF=/data/ban-m/a_thaliana/genome/Arabidopsis_thaliana.TAIR10.42.gff3


## ===== Output =====
FASTA=./data/plant_mitochondria_from_genbank.fa
TAB=./data/plant_mitochondria_from_genbank.tab


## ==== Main proc ====

python3 ./script/extract_plant_mitochondria_from_genbank.py ${GENBANK} ${FASTA} ${TAB}
cat ${ARABIDOPSIS} >> ${FASTA}
cargo run --release --bin convert_gff ${ARABIDOPSIS_GFF} >> ${TAB}
