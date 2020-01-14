#!/bin/bash
Rscript ./scripts/plot_length_and_accuracy.R \
         ./result/length_and_accuracy_default.tsv length_and_accuracy_default
Rscript ./scripts/plot_length_and_accuracy.R \
         ./result/length_and_accuracy_pacbio.tsv length_and_accuracy_pacbio

# Rscript ./scripts/plot_few_sample.R \
#         ./result/few_sample_default.tsv few_sample_default
# Rscript ./scripts/plot_few_sample.R \
#         ./result/few_sample_pacbio.tsv few_sample_pacbio

Rscript ./scripts/plot_neutral_prediction.R
Rscript ./scripts/plot_coverage_and_likelihood.R
Rscript ./scripts/plot_coverage_and_contained.R
Rscript ./scripts/plot_number_of_edge_and_node.R

