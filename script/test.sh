#!/bin/bash
set -ue
OUTPUT=${PWD}/test
mkdir -p ${OUTPUT}
# LEN=100000
# cargo run --release --bin create_mock_genomes -- ${LEN} 0.001 2 13223 \
#       ${OUTPUT}/testseq
# badread simulate \
#         --reference ${OUTPUT}/testseq_contigs.fa \
#         --quantity 50x --error_model pacbio \
#         --seed 10\
#         --qscore_model pacbio --identity 90,95,3 \
#         --junk_reads 0 --random_reads 0 --chimeras 0 \
#         --length 15000,5000 > ${OUTPUT}/testseq_reads.fq
# cat ${OUTPUT}/testseq_reads.fq |\
#     paste - - - - | cut -f 1,2 |\
#     sed -e 's/@/>/g' | tr '\t' '\n' \
#                           > ${OUTPUT}/testseq_reads.fa
# READ=${OUTPUT}/testseq_reads.fa
# REFERENCE=${OUTPUT}/testseq_reference.fa
# qsub -o ./logfiles/test.log -j y -S /bin/bash -cwd -pe smp 24 -V \
#      ./script/disassembly.sh \
#      ${REFERENCE} \
#      ${READ} \
#      ${PWD}/test/test 2 1800 24

# OUTPUT=${PWD}/test/test
# READ=${PWD}/test/testseq_reads.fa
# REFERENCE=${PWD}/test/testseq_reference.fa
# MIN_CLUSTER=2
# CORES=12
# LIMIT=1800
# ${PWD}/target/release/mmmm \
#       decompose --alignments ${OUTPUT}/last_db/alignments.tab --output ${OUTPUT} \
#       --reads ${READ} --contigs ${REFERENCE} \
#       --self_alignments ${OUTPUT}/last_db/self.tab \
#       --cluster_num ${MIN_CLUSTER} --threads ${CORES} \
#       --limit ${LIMIT}\
#       --resume ${PWD}/test/198.json \
#       -vv

OUTPUT=${PWD}/result/synthetic_data/150_001
READ=${PWD}/data/synthetic_data/150_001/150_001_reads.fa
REFERENCE=${PWD}/data/synthetic_data/150_001/150_001_reference.fa
MIN_CLUSTER=3
CORES=12
LIMIT=1213
cargo build --release 
${PWD}/target/release/mmmm \
      decompose --alignments ${OUTPUT}/last_db/alignments.tab --output ${OUTPUT} \
      --reads ${READ} --contigs ${REFERENCE} \
      --self_alignments ${OUTPUT}/last_db/self.tab \
      --cluster_num ${MIN_CLUSTER} --threads ${CORES} \
      --limit ${LIMIT}\
      --resume ${OUTPUT}/encoded_reads.json\
      -vv 2> 150_001.log

OUTPUT=${PWD}/result/synthetic_data/100_01
READ=${PWD}/data/synthetic_data/100_01/100_01_reads.fa
REFERENCE=${PWD}/data/synthetic_data/100_01/100_01_reference.fa
MIN_CLUSTER=3
CORES=12
LIMIT=1213
cargo build --release 
${PWD}/target/release/mmmm \
      decompose --alignments ${OUTPUT}/last_db/alignments.tab --output ${OUTPUT} \
      --reads ${READ} --contigs ${REFERENCE} \
      --self_alignments ${OUTPUT}/last_db/self.tab \
      --cluster_num ${MIN_CLUSTER} --threads ${CORES} \
      --limit ${LIMIT}\
      --resume ${OUTPUT}/encoded_reads.json\
      -vv 2> 100_01.log



