#!/bin/bash
set -ue
threads=6
for i in `seq 0 20 `
do
    OUTPUT=${PWD}/benchmark/benchmark_${i}.txt
    OUTLOG_PREFIX=${PWD}/benchmark/last_decompose_num_poa_${i}
    mkdir -p ${PWD}/benchmark
    COVERAGE=0
    PROBS="0.5 0.5"
    seed=1213
    seed=$(( $seed * $i ))
    ${PWD}/target/release/test_last_decompose_multiple_varying_errorate \
          ${seed} ${threads} ${PROBS} > ${OUTPUT} \
          2> ${OUTLOG_PREFIX}.log
done

cat ${PWD}/result/benchmark/benchmark*.txt |\
    grep RESULT > ${PWD}/result/benchmark/benchmark.tsv
