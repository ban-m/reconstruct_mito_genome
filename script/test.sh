#!/bin/bash
COVERAGE=0
PROBS="0.5 0.5"
seed=1213
threads=24
seed=$(( $seed * 2 ))
RUST_LOG=debug ${PWD}/target/release/test_last_decompose_multiple_varying_errorate \
        ${seed} ${threads} ${PROBS} > temp.tsv \
        2> test.log
cat test | grep RESULT >> out.tsv
