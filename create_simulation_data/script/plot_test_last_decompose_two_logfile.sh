#!/bin/bash
for frac in skew even
do
    for readnum in 100 300
    do
        cat ./logfiles/last_decompose_vb_prior.${frac}.${readnum}.log |\
            grep Summary |\
            grep -v ID > ./result/last_decompose_vb_prior_${frac}_${readnum}_log.tsv
        Rscript --vanilla --slave ./script/plot_test_last_decompose_two_logfile.R \
                ./result/last_decompose_vb_prior_${frac}_${readnum}_log.tsv\
                last_decompose_vb_prior_${frac}_${readnum}
    done
done

for frac in skew even
do
    for readnum in 100 300
    do
        cat ./logfiles/last_decompose_multi_vb_prior.${frac}.${readnum}.log |\
            grep Summary |\
            grep -v ID > ./result/last_decompose_multi_vb_prior_${frac}_${readnum}_log.tsv
        Rscript --vanilla --slave ./script/plot_test_last_decompose_two_logfile.R \
                ./result/last_decompose_multi_vb_prior_${frac}_${readnum}_log.tsv\
                last_decompose_multi_vb_prior_${frac}_${readnum}
    done
done
