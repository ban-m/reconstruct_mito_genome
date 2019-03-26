#$ -S /bin/bash
#$ -N debug
#$ -cwd
#$ -pe smp 12
#$ -e ./logfiles/start_and_stop.log
#$ -o ./logfiles/start_and_stop.out
#$ -V
#$ -m e
set -ue
function procedure() {
    BAM=$1
    OUTPUT=$2
    REFERENCE=$3
    READ=$4
    if [ ! -f ${1} ]
    then
        echo "ERROR" ${1} " does not exist"
        return 1
    fi
    
    if [ ! -f ${3} ]
    then
        echo "ERROR" ${3} " does not exist"
        return 1
    fi

    if [ ! -f ${4} ]
    then
        echo "ERROR" ${4} " does not exist"
        return 1
    fi

    cargo run --release \
          --bin start_and_stopping_read \
          -- ${BAM} \
          > ${OUTPUT}
    cargo run --release \
          --bin peak_call_from_start_stop \
          -- ${OUTPUT} \
          > ${OUTPUT%.tsv}.peaks.tsv
    Rscript --vanilla --slave ./script/start_and_stop_read_plot.R ${OUTPUT} ${OUTPUT%.tsv}.peaks.tsv
}

export -f procedure
SEQUEL_ROOT=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback
SEQUEL_READ=/data/ban-m/a_thaliana/sequel_reads/sequel1_filter_dedupe.fq
procedure ${SEQUEL_ROOT}/guided_asm_sequel_canu.contigs.fasta.mapback.bam \
          ${SEQUEL_ROOT}/guided_asm_sequel_canu.contigs.fasta.mapback.tsv \
          /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/canu/guided_asm_sequel_canu.contigs.fasta \
          ${SEQUEL_READ}

procedure ${SEQUEL_ROOT}/scaffolds.fasta.mapback.bam \
          ${SEQUEL_ROOT}/scaffolds.fasta.mapback.tsv \
          /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/flye/scaffolds.fasta \
          ${SEQUEL_READ}

procedure ${SEQUEL_ROOT}/wtdbg2_sequel.ctg.1.fa.mapback.bam \
          ${SEQUEL_ROOT}/wtdbg2_sequel.ctg.1.fa.mapback.tsv \
          /grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/wtdbg2/wtdbg2_sequel.ctg.1.fa \
          ${SEQUEL_READ}



# ONT_READ=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/filtered_reads.fq
# ONT_ROOT=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/mapback
# find ./result/ -name "*nanopore*.fa*" |\
#     parallel procedure ${ONT_ROOT}/{/}.mapback.bam \
#              ${ONT_ROOT}/{/}.mapback.tsv \
#              {} \
#              ${ONT_READ}

