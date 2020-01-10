#!/bin/bash
# qsub ${PWD}/script/yellowgate.job \
#      ${PWD}/data/mock_genome_read.fa\
#      ${PWD}/data/mock_genome_ref.fa \
#      ${PWD}/result/mock_genome_yg

# qsub ./script/yellowgate.job \
#      ${PWD}/data/easy/mock_genome_read.fa \
#      ${PWD}/data/easy/mock_genome_ref.fa \
#      ${PWD}/result/short_easy/


# qsub ./script/yellowgate.job \
#      ${PWD}/data/short_hard/reads.fa \
#      ${PWD}/data/short_hard/mock_genome_ref.fa \
#      ${PWD}/result/short_hard

# qsub ./script/yellowgate.job \
#      ${PWD}/data/middle_hard/reads.fa \
#      ${PWD}/data/middle_hard/mock_genome_ref.fa \
#      ${PWD}/result/middle_hard

# qsub ./script/yellowgate.job \
#      ${PWD}/data/long_hard/reads.fa\
#      ${PWD}/data/long_hard/mock_genome_ref.fa\
#      ${PWD}/result/long_hard

# qsub ./script/yellowgate.job \
#      ${PWD}/data/short_extreme/reads.fa \
#      ${PWD}/data/short_extreme/mock_genome_ref.fa \
#      ${PWD}/result/short_extreme

qsub ./script/yellowgate.job \
     ${PWD}/data/middle_extreme/reads.fa \
     ${PWD}/data/middle_extreme/mock_genome_ref.fa \
     ${PWD}/result/middle_extreme

# qsub ./script/yellowgate.job \
#      ${PWD}/data/long_extreme/reads.fa \
#      ${PWD}/data/long_extreme/mock_genome_ref.fa \
#      ${PWD}/result/long_extreme

# qsub ./script/yellowgate.job\
#      ${PWD}/data/complex/read_complex1.fa\
#      ${PWD}/data/complex/reference.fa \
#      ${PWD}/result/complex1

# qsub ./script/yellowgate.job\
#      ${PWD}/data/complex/read_complex2.fa\
#      ${PWD}/data/complex/reference.fa \
#      ${PWD}/result/complex2

