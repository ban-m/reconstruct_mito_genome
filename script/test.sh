#!/bin/bash
ROOT=https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly
for accession in pacbio ler col0_1106_exp2 cvi an1 c24 kyoto sha eri
do
    reads=${ROOT}/filtered_reads/${accession}.fasta
    result=${ROOT}/${accession}.tar.gz
    circos=${ROOT}/viewer/${accession}/circos.html
    linear=${ROOT}/viewer/${accession}/linear.html
    no_merge=${ROOT}/viewer/${accession}/no_merge.html
    echo -e "|${accession}|[Reads](${reads})|[Result](${result})|[Circos](${circos})|[Liner](${linear})|[NoMerge](${no_merge})|"
done
