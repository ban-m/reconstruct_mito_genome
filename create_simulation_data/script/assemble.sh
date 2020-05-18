#!/bin/bash

qsub -j y -o ${PWD}/logfiles/complex1.log ./script/test.job \
        ${PWD}/data/complex/read_complex1.fa \
        ${PWD}/result/complex1

qsub -j y -o ${PWD}/logfiles/complex2.log ./script/test.job \
        ${PWD}/data/complex/read_complex2.fa \
        ${PWD}/result/complex2

