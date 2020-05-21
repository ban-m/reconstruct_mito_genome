#!/bin/bash
cat ${PWD}/logfiles/two_test_gibbs.log | rg FEATURE > ${PWD}/result/lks_two.tsv
cat ${PWD}/logfiles/multi_test.log | rg FEATURE > ${PWD}/result/lks_six.tsv
Rscript --vanilla --slave ${PWD}/script/plot_clustering.R 
