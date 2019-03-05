#!/bin/bash



mkdir -p /grid/ban-m/arabidopsis_thaliana/nanopore/coverage_asm/
qsub ./script/read_coverage_assembly.job \
     /data/ban-m/a_thaliana/ONT/ERR2173373.fastq \
     /grid/ban-m/arabidopsis_thaliana/nanopore/coverage_asm \
     readcov_asm_nanopore_canu \
     ONT

# mkdir -p /grid/ban-m/arabidopsis_thaliana/sequel/coverage_asm/
# qsub ./script/read_coverage_assembly.job \
#      /data/ban-m/a_thaliana/sequel_reads/sequel1_filter_dedupe.fq \
#      /grid/ban-m/arabidopsis_thaliana/sequel/coverage_asm \
#      readcov_asm_sequel_canu \
#      PacBio

