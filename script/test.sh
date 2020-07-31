OUTPUT=${PWD}/test
mkdir -p ${OUTPUT}
# TESTNUM=200
# COVERAGE=0
# SEED=324
# DIV_RATE=0.001
# PROBS="0.5 0.5"
# for seed in 12 348309 49232 23923 432
# do
#     RUST_LOG=debug cargo run --release --bin test_last_decompose_multiple -- \
#             ${TESTNUM} ${COVERAGE} \
#             ${seed} ${DIV_RATE} ${PROBS}  \
#             >> ${OUTPUT}/test.out 2>> ${OUTPUT}/test.log
# done
# RUST_LOG=debug ${PWD}/target/release/test_last_decompose_multiple_varying_errorate \
#         ${SEED} ${PROBS} > ${OUTPUT}/varying.out 2> ${OUTPUT}/varying.log

LEN=100000
cargo run --release --bin create_mock_genomes -- ${LEN} 0.001 2 13223 \
      ${OUTPUT}/testseq
badread simulate \
        --reference ${OUTPUT}/testseq_contigs.fa \
        --quantity 50x --error_model pacbio \
        --seed 10\
        --qscore_model pacbio --identity 90,95,3 \
        --junk_reads 0 --random_reads 0 --chimeras 0 \
        --length 15000,5000 > ${OUTPUT}/testseq_reads.fq
cat ${OUTPUT}/testseq_reads.fq |\
    paste - - - - | cut -f 1,2 |\
    sed -e 's/@/>/g' | tr '\t' '\n' \
                          > ${OUTPUT}/testseq_reads.fa
READ=${OUTPUT}/testseq_reads.fa
REFERENCE=${OUTPUT}/testseq_reference.fa
qsub -o ./logfiles/test.log -j y -S /bin/bash -cwd -pe smp 24 -V \
     ./script/disassembly.sh \
     ${REFERENCE} \
     ${READ} \
     ${PWD}/test/test 2 1800 24

