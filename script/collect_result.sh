#!/bin/bash
# set -ue 
mkdir -p ${PWD}/result/viewer
cp ${PWD}/script/circos.js ${PWD}/result/viewer/circos.js
cp ${PWD}/result/pacbio/viewer/style.css ${PWD}/result/viewer/style.css
cp ${PWD}/script/linear.js ${PWD}/result/viewer/linear.js
ROOT="./.."
for accession in pacbio an1 c24 cvi eri kyo ler sha col0_1106_exp2 
do
    DIST=${PWD}/result/viewer/${accession}
    mkdir -p ${PWD}/result/viewer/${accession}
    # Circos plot
    cat ${PWD}/result/${accession}/viewer/circos.html |\
        sed "s+/viewer/style.css+${ROOT}/style.css+g" |\
        sed "s+/viewer/circos.js+${ROOT}/circos.js+g" |\
        sed "s+/viewer/data.json+./data_circos.json+g" |\
        sed "s+/viewer/repeats.json+./repeats.json+g"> ${DIST}/circos.html
    cp ${PWD}/result/${accession}/viewer/data.json ${DIST}/data_circos.json
    cp ${PWD}/result/${accession}/viewer/repeats.json ${DIST}/repeats.json
    # Circos plot no_merge.
    cat ${PWD}/result/${accession}/viewer/circos.html |\
        sed "s+/viewer/style.css+${ROOT}/style.css+g" |\
        sed "s+/viewer/circos.js+${ROOT}/circos.js+g" |\
        sed "s+/viewer/data.json+./data_no_merge.json+g" |\
        sed "s+/viewer/repeats.json+./repeats.json+g"> ${DIST}/no_merge.html
    cp ${PWD}/result/${accession}/no_merge/viewer/data.json ${DIST}/data_no_merge.json
    cat ${PWD}/result/${accession}/viewer/linear.html |\
        sed "s+/viewer/style.css+${ROOT}/style.css+g" |\
        sed "s+/viewer/linear.js+${ROOT}/linear.js+g" |\
        sed "s+/viewer/read_data.json+./read_data.json+g" |\
        sed "s+/viewer/contig_alns.json+./contig_alns.json+g" > ${DIST}/linear.html
    cp ${PWD}/result/${accession}/viewer/read_data.json ${DIST}/read_data.json
    cp ${PWD}/result/${accession}/viewer/contig_alns.json ${DIST}/contig_alns.json
done
