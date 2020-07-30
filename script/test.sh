OUTPUT=${PWD}/data/test
mkdir -p ${OUTPUT}
LEN=100000
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

READ=${OUTPUT}/testseq_reads.fa
REFERENCE=${OUTPUT}/testseq_reference.fa
qsub -o ./logfiles/test.log -j y -S /bin/bash -cwd -pe smp 24 -V \
     ./script/disassembly.sh \
     ${REFERENCE} \
     ${READ} \
     ${PWD}/test/test 2 40000 24

