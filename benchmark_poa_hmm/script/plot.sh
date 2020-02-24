#!/bin/bash
# Rscript --slave --vanilla \
#         ./script/plot_few_sample.R ./result/few_sample_default.tsv few_sample_default
# Rscript --slave --vanilla \
#         ./script/plot_few_sample.R ./result/few_sample_pacbio.tsv few_sample_pacbio

# Rscript --slave --vanilla \
#         ./script/plot_length_and_accuracy.R \
#         ./result/length_and_accuracy_pacbio.tsv length_and_accuracy_pacbio

Rscript --slave --vanilla \
        ./script/plot_length_and_accuracy.R \
        ./result/length_and_accuracy_default.tsv length_and_accuracy_default
