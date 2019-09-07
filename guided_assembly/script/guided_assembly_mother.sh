#!/bin/bash
#$ -S /bin/bash
#$ -N MitoAsm
#$ -cwd
#$ -pe smp 24
#$ -V
#$ -e ./logfiles/mitoasm.log
#$ -o ./logfiles/mitoasm.out
#$ -m e
# -M banmasutani@gmail.com ## To send a mail when ended
# -t 1:n ## For array job
set -ue

### Variables

ONT=/data/ban-m/a_thaliana/ONT/ERR2173373.fastq
SEQUEL=/data/ban-m/a_thaliana/sequel_reads/sequel1_filter_dedupe.fq
REFERENCE=/grid/ban-m/arabidopsis_thaliana/genome/GCF_000001735.4_TAIR10.1_genomic.fna

ONT_OUTPUT=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/
SEQUEL_OUTPUT=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/

ONT_FILTERED=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/filtered_reads.fq
SEQUEL_FILTERED=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fq


### ==== Preprocess =====

rm -rf ${ONT_OUTPUT} ${SEQUEL_OUTPUT}
mkdir -p ${ONT_OUTPUT}
mkdir -p ${SEQUEL_OUTPUT}

#### ==== Filter ONT read =====

## Note that in RefSeq 2019, the mitochondrion sequence is called NC_037304.1.
TEMP=${RANDOM}.sam
minialign -P -t 12 -x ont.r9 ${REFERENCE} ${ONT} | \
    awk -F'\t' '($5 >= 50 && $3 == "NC_037304.1"){print $0}'\
        > ${TEMP}
cat ${TEMP} | cut -f 1 | cargo run --release \
                               --bin main \
                               -- ${ONT} > ${ONT_FILTERED}

#### ==== Filter Sequel Read =====
minialign -P -t 12 -x pacbio ${REFERENCE} ${SEQUEL} | \
    awk -F'\t' '($5 >= 50 && $3 == "NC_037304.1"){print $0}' \
        >  ${TEMP}
cat ${TEMP} | cut -f 1 | cargo run --release \
                               --bin main \
                               -- ${SEQUEL} > ${TEMP}.fq
cargo run --release --bin clip_self_chimera -- ${TEMP}.fq > ${SEQUEL_FILTERED}
rm ${TEMP} ${TEMP}.fq

#### ==== Calculate Stats for each readset =====

# echo -e "Name\t# of read\tmax length\tmin length\tN50\tTotal\tMean length\tEst.cov" > ./result/whole_read_summary.tsv
# read_statistics ${ONT} 400000 ONT_whole >> ./result/whole_read_summary.tsv
# read_statistics ${SEQUEL} 400000 sequel_whole >> ./result/whole_read_summary.tsv

# echo -e "Name\t# of read\tmax length\tmin length\tN50\tTotal\tMean length\tEst.cov" > ./result/read_summary.tsv
# read_statistics ${ONT_FILTERED} 400000 ONT >> ./result/read_summary.tsv
# read_statistics ${SEQUEL_FILTERED} 400000 sequel >> ./result/read_summary.tsv

# K=32
# jellyfish count -m 32 -s 100M -t 12 -o ./result/ont_${K}mer_freq.jf --bf-size 10G ${ONT_FILTERED}
# jellyfish histo -o ./result/ont_${K}mer_freq.tsv -t 12 ./result/ont_${K}mer_freq.jf
# jellyfish count -m 32 -s 100M -t 12 -o ./result/sequel_${K}mer_freq.jf --bf-size 10G ${SEQUEL_FILTERED}
# jellyfish histo -o ./result/sequel_${K}mer_freq.tsv -t 12 ./result/sequel_${K}mer_freq.jf

#### ==== Canu =====
rm -rf ${ONT_OUTPUT}/canu
mkdir -p ${ONT_OUTPUT}/canu
qsub -o ./logfiles/canu_nanopore.out -e ./logfiles/canu_nanopore.log \
     ./script/canu.job \
     ${ONT_OUTPUT}/canu \
     guided_asm_nanopore_canu \
     -nanopore-raw \
     ${ONT_FILTERED} 

rm -rf ${SEQUEL_OUTPUT}/canu/
mkdir -p ${SEQUEL_OUTPUT}/canu/
qsub -o ./logfiles/canu_sequel.out -e ./logfiles/canu_sequel.log \
     ./script/canu.job \
     ${SEQUEL_OUTPUT}/canu \
     guided_asm_sequel_canu \
     -pacbio-raw \
     ${SEQUEL_FILTERED}

#### ==== Ra ======
rm -rf ${ONT_OUTPUT}/ra
mkdir -p ${ONT_OUTPUT}/ra
qsub -o ./logfiles/ra_nanopore.out -e ./logfiles/ra_nanopore.log \
     ./script/ra.job \
     ${ONT_OUTPUT}/ra \
     ont \
     ${ONT_FILTERED} \
     guided_asm_nanopore_ra.fa

rm -rf ${SEQUEL_OUTPUT}/ra/
mkdir -p ${SEQUEL_OUTPUT}/ra/
qsub -o ./logfiles/ra_sequel.out -e ./logfiles/ra_sequel.log \
     ./script/ra.job \
     ${SEQUEL_OUTPUT}/ra \
     pb \
     ${SEQUEL_FILTERED} \
     guided_asm_sequel_ra.fa


#### ==== Wtdbg2 =====
rm -rf ${ONT_OUTPUT}/wtdbg2
mkdir -p ${ONT_OUTPUT}/wtdbg2
qsub -o ./logfiles/wtdbg2_nanopore.out -e ./logfiles/wtdbg2_nanopore.log \
     ./script/wtdbg2.job \
     ${ONT_OUTPUT}/wtdbg2 \
     ont \
     ${ONT_FILTERED} \
     wtdbg2_nanopore \
     map-ont 

rm -rf ${SEQUEL_OUTPUT}/wtdbg2
mkdir -p ${SEQUEL_OUTPUT}/wtdbg2
qsub -o ./logfiles/wtdbg2_sequel.out -e ./logfiles/wtdbg2_sequel.log \
     ./script/wtdbg2.job \
     ${SEQUEL_OUTPUT}/wtdbg2 \
     sequel \
     ${SEQUEL_FILTERED} \
     wtdbg2_sequel \
     map-pb 

#### ===== Flye =======
rm -rf ${ONT_OUTPUT}/flye
mkdir -p ${ONT_OUTPUT}/flye
qsub -o ./logfiles/flye_nanopore.out -e ./logfiles/flye_nanopore.log \
     ./script/flye.job \
     nano-raw \
     ${ONT_FILTERED} \
     ${ONT_OUTPUT}/flye \
     
rm -rf ${SEQUEL_OUTPUT}/flye
mkdir -p ${SEQUEL_OUTPUT}/flye
qsub -o ./logfiles/flye_sequel.out -e ./logfiles/flye_sequel.log \
     ./script/flye.job \
     pacbio-raw \
     ${SEQUEL_FILTERED} \
     ${SEQUEL_OUTPUT}/flye


##### ===== FALCON =======


## Postprocess

bash ./script/guided_assembly_postprocess.sh
