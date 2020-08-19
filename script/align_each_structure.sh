#!/bin/bash
set -ue
ROOT=${PWD}
OUTPUT=${PWD}/result/alignments/
OUTFILE=${OUTPUT}/summary.tsv
DATA=${PWD}/result
mkdir -p ${OUTPUT}
cd ${OUTPUT}
CORES=12
echo -e "Target\tQuery\tCoverRate" > ${OUTFILE}

for reference in an1 c24 cvi eri kyo ler sha col0_1106_exp2
do
    for query in an1 c24 cvi eri kyo ler sha col0_1106_exp2
    do
        if [ ${reference} = ${query} ]
        then
            echo "${reference} and ${query} are the same."
            echo -e "${reference}\t${query}\t1" >> ${OUTFILE}
        else
            ref_seq=${DATA}/${reference}/multipartite.fasta
            query_seq=${DATA}/${query}/multipartite.fasta
            lastdb -R00 -Q0 ${reference} ${ref_seq} 
            last-train -Q0 ${reference} ${query_seq} > score.matrix
            lastal -P${CORES} -f tab -R00 -Q0 -p score.matrix \
                   ${reference} ${query_seq} \
                   > ${query}_to_${reference}.tab
            cover_rate=$(cargo run --release --bin overlaps -- ${query}_to_${reference}.tab ${ref_seq})
            echo -e "${reference}\t${query}\t${cover_rate}" >> ${OUTFILE}
        fi
    done
done

sed -i -e 's/_1106_exp2//g' ${OUTFILE} 
cd ${ROOT}
Rscript --vanilla --slave ./script/plot_alignment.R ${OUTFILE} aln_heatmap
