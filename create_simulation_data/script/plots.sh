#!/bin/bash
cat ${PWD}/result/last_decompose_num_poa_* | rg RESULT > ${PWD}/result/last_decompose_num_poa.tsv
Rscript --vanilla --slave ${PWD}/script/plot_errors.R ${PWD}/result/last_decompose_num_poa.tsv last_decompose_num_poa

