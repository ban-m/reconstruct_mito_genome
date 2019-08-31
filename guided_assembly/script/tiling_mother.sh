#!/bin/bash
qsub ./script/tiling_by_stiff_sequence.job \
     /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fq \
     /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/canu/guided_asm_sequel_canu.contigs.fasta \
     /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback/guided_asm_sequel_canu.contigs.fasta.mapback.peaks.tsv \
     /grid/ban-m/arabidopsis_thaliana/sequel/tiling \
     canu
