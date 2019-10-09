#!/bin/bash
qsub ./scripts/workflow.job \
     /data/ban-m/a_thaliana/sequel_reads/ERR341582.fasta /grid/ban-m/arabidopsis_thaliana/sequel/assemble/mpi \
     /grid/ban-m/arabidopsis_thaliana/genome/GCF_000001735.4_TAIR10.1_genomic.fna

     
