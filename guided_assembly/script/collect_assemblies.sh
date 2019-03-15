#!/bin/bash
mkdir -p result/assemblies/

cp /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/canu/guided_asm_sequel_canu.contigs.fasta \
   ./result/assemblies/sequel_canu.fasta
cp /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/ra/guided_asm_sequel_ra.fa \
   ./result/assemblies/sequel_ra.fasta
cp /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/wtdbg2/wtdbg2_sequel.ctg.1.fa \
   ./result/assemblies/sequel_wtdbg2.fasta
cp /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/flye/scaffolds.fasta \
   ./result/assemblies/sequel_flye.fasta


cp /grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/canu/guided_asm_nanopore_canu.contigs.fasta \
   ./result/assemblies/nanopore_canu.fasta
cp /grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/ra/guided_asm_nanopore_ra.fa \
   ./result/assemblies/nanopore_ra.fasta
cp /grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/wtdbg2/wtdbg2_nanopore.ctg.1.fa \
   ./result/assemblies/nanopore_wtdbg2.fasta
cp /grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/flye/scaffolds.fasta \
   ./result/assemblies/nanopore_flye.fasta
