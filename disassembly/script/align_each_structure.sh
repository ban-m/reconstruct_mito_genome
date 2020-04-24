#!/bin/bash
ROOT=${PWD}
OUTPUT=${PWD}/result/alignments/
REFERENCE=${PWD}/result/pacbio/multipartite.fasta
DATA=${PWD}/result
mkdir -p ${OUTPUT}
cd ${OUTPUT}
lastdb -R00 -Q0 reference ${REFERENCE}
CORES=12

# for accession in an1 c24 cvi eri kyo ler sha col0_1106_exp2 tal6111_1106_exp2 tal6144_1115_exp2 tal61up63_1106_exp2 tal6226_1115_exp2
# do
#     QUERY=${DATA}/${accession}/multipartite.fasta
#     last-train -Q0 reference ${QUERY} > score.matrix
#     lastal -P${CORES} -f tab -R00 -Q0 -p score.matrix reference ${QUERY} |\
#         rg -v "#" |\
#         awk -v accession=${accession} 'BEGIN{OFS="\t"}{$7 = accession ; print $0}' \
#             > ${accession}.tab
#     cargo run --release --bin overlaps -- ${REFERENCE} ${accession}.tab |\
#         awk -v accession=${accession} 'BEGIN{OFS="\t"}{print $0, accession}' \
#             >> occurrences.tsv
# done
cd ${ROOT}
Rscript --vanilla --slave ./script/plot_alignment.R ${OUTPUT}/occurrences.tsv aln_heatmap
