#!/bin/bash
set -ue

if [ $# -ne 5 ]
then
    echo "Error: the length of the input is not 5. please give the arguments correctly to the program." 1>&2
    exit 1
fi


#### ----- Sanity check ------

INPUT=$1
OUTPUT_DIR=$2
PREFIX=$3
DATABASE=$4
ANNOTATION_FILE=$5

if [ ! -e ${INPUT} ]
then
    echo "Error: the input ${INPUT} does not exists." 1>&2
    exit 1
fi

mkdir -p $OUTPUT_DIR

if [ ! -d ${OUTPUT_DIR} ]
then
    echo "Error: can not create the ${OUTPUT} dir." 1>&2
    exit 1
fi

if [ ! -e ${DATABASE} ]
then
    echo "Error: could not find ${DATABASE}." 1>&2
    exit 1
fi

if [ ! -e ${ANNOTATION_FILE} ]
then
    echo "Error: could not find ${ANNOATION_FILE}" 1>&2
    exit 1
fi


### ----- LAST alignments -------
# echo "lastdb start..."
# mkdir -p ${OUTPUT_DIR}/last/
# lastdb -P 12 -R11 -uYASS -Q0 ${OUTPUT_DIR}/last/lastdb_${PREFIX} ${DATABASE}
# echo "lastal start..."
# lastal -f TAB ${OUTPUT_DIR}/last/lastdb_${PREFIX} ${INPUT} > ${OUTPUT_DIR}/last/${PREFIX}_todb.tab
# echo  "Success."


### ----- tRNAscan-SE ------
# echo "tRNAscan-SE version2.0..."
# rm -r ${OUTPUT_DIR}/tRNAscan-SE/
# mkdir -p ${OUTPUT_DIR}/tRNAscan-SE/
# tRNAscan-SE -O ${INPUT} -o ${OUTPUT_DIR}/tRNAscan-SE/${PREFIX}_trnascan.txt \
#             -f ${OUTPUT_DIR}/tRNAscan-SE/${PREFIX}_trnascan.structure.dat \
#             -m ${OUTPUT_DIR}/tRNAscan-SE/${PREFIX}_trnascan.summary.txt 
# echo "Success."

### ---- Summarize ------
echo "Convert last TAB file into Annotation JSON..."
cargo run --release --bin lasttab_to_json --  \
      ${OUTPUT_DIR}/last/${PREFIX}_todb.tab ${ANNOTATION_FILE} \
      > ${OUTPUT_DIR}/${PREFIX}_last.json
# echo "Success."
# echo "Convert tRNAscanSE file into Annotation JSON..."
# cargo run --release --bin trnascan_to_json -- \
#       ${OUTPUT_DIR}/tRNAscan-SE/${PREFIX}_trnascan.txt \
#       > ${OUTPUT_DIR}/${PREFIX}_trnascan.json
# echo "Success."

### ---- Visualize -----
# echo "Genetating HTML/SGV figures..."
# cargo run --release --bin convert_fasta -- ${INPUT} > ${OUTPUT_DIR}/${PREFIX}_contig.json
# cat ./data/template.html |\
#     sed -e "s+LAST_JSON+${OUTPUT_DIR}/${PREFIX}_last.json+g" |\
#     sed -e "s+tRNASE_JSON+${OUTPUT_DIR}/${PREFIX}_trnascan.json+g" |\
#     sed -e "s+FASTA_JSON+${OUTPUT_DIR}/${PREFIX}_contig.json+g" |\
#     > ${OUTPUT_DIR}/viewer.html
# cp ./data/drawer.js ${OUTPUT_DIR}/
# cp ./data/style.css ${OUTPUT_DIR}/
# echo "Success. Opne ${OUTPUT_DIR}/viewer.html in a modern web browser."
