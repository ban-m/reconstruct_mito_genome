#!/bin/bash
#$ -S /bin/bash
#$ -N PostProcess
#$ -cwd
#$ -pe smp 24
#$ -e ./logfiles/guided_assembly_postprocess.log
#$ -o ./logfiles/guided_assembly_postprocess.out
#$ -V
#$ -m e

set -ue

### ==== Result =====
ONT=/data/ban-m/a_thaliana/ONT/ERR2173373.fastq
ONT_FILTERED=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/filtered_reads.fq

ONT_RA=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/ra/guided_asm_nanopore_ra.fa
ONT_CANU=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/canu/guided_asm_nanopore_canu.contigs.fasta
ONT_WTDBG2=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/wtdbg2/wtdbg2_nanopore.ctg.1.fa
ONT_RA_GFA=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/ra/rala_assembly_graph.gfa
ONT_CANU_GFA=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/canu/guided_asm_nanopore_canu.contigs.gfa

SEQUEL=/data/ban-m/a_thaliana/sequel_reads/sequel1_filter_dedupe.fq
SEQUEL_FILTERED=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/filtered_reads.fq

SEQUEL_RA=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/ra/guided_asm_sequel_ra.fa
SEQUEL_CANU=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/canu/guided_asm_sequel_canu.contigs.fasta
SEQUEL_WTDBG2=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/wtdbg2/wtdbg2_sequel.ctg.1.fa
SEQUEL_FLYE=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/flye/scaffolds.fasta
SEQUEL_RA_GFA=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/ra/rala_assembly_graph.gfa
SEQUEL_CANU_GFA=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/canu/guided_asm_sequel_canu.contigs.gfa
SEQUEL_FLYE_GFA=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/flye/assembly_graph.gfa

### ==== Post process of guided assembly =====

#### ==== Mapping back =====

ONT_OUTPUT=/grid/ban-m/arabidopsis_thaliana/nanopore/guided_asm/mapback
mkdir -p ${ONT_OUTPUT}
for contigs in ${ONT_RA} ${ONT_CANU} ${ONT_WTDBG2}
do
    if [ -e  ${contigs} ]
    then
        prefix=${ONT_OUTPUT}/${contigs##*/}
        minimap2 -t 24 -a -x map-ont ${contigs} ${ONT_FILTERED} | \
            samtools view -bhS | \
            samtools sort -@ 12 -m 5G -O BAM > ${prefix}.mapback.raw.bam
        cargo run --release --bin filtering_out ${prefix}.mapback.raw.bam ${prefix}.mapback.bam
        samtools stats ${prefix}.mapback.bam > ${prefix}.mapback.stats
        samtools index ${prefix}.mapback.bam
    fi
done

SEQUEL_OUTPUT=/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback
mkdir -p ${SEQUEL_OUTPUT}
for contigs in ${SEQUEL_RA} ${SEQUEL_WTDBG2} ${SEQUEL_CANU} ${SEQUEL_FLYE}
do
    prefix=${SEQUEL_OUTPUT}/${contigs##*/}
    minimap2 -t 24 -a -x map-pb ${contigs} ${SEQUEL_FILTERED} | \
        samtools view -hbS | \
        samtools sort -@12 -m 5G -O BAM > ${prefix}.mapback.raw.bam
    cargo run --release --bin filtering_out ${prefix}.mapback.raw.bam ${prefix}.mapback.bam
    samtools stats ${prefix}.mapback.bam > ${prefix}.mapback.stats
    samtools index ${prefix}.mapback.bam
done

# ### ===== Collecting results =======

# mkdir -p result
# cp $ONT_RA ./result/
# cp $ONT_RA_GFA ./result/guided_asm_ra_nanopore.gfa
# cp $ONT_WTDBG2 ./result/
# cp $ONT_CANU ./result/
# cp $ONT_CANU_GFA ./result/

# cp $SEQUEL_RA ./result/
# cp $SEQUEL_RA_GFA ./result/guided_asm_ra_sequel.gfa
# cp $SEQUEL_WTDBG2 ./result/
# cp $SEQUEL_CANU ./result/
# cp $SEQUEL_CANU_GFA ./result/


# ### ==== Tar gzipping ====
# cd result
# tar -czvf guided_assemblies.tar.gz *.fa *.fasta *.gfa
# tar -czvf mapback.tar.gz ${ONT_OUTPUT}/* ${SEQUEL_OUTPUT}/* 
# cd ../

# ### ==== Summarize =====
# rm -f ./result/assembly_stats.dat ./result/assembly_stats_mini.tsv
# echo -e "Name\tN50\tMean\t# of contig\t# of N\t$ of Gaps" > ./result/assembly_stats_mini.tsv
# for file in ./result/*.fa
# do
#     echo "======" $file "======" >> ./result/assembly_stats.dat
#     assembly-stats $file >> ./result/assembly_stats.dat
#     echo "" >> ./result/assembly_stats.dat
#     assembly-stat $file ${file##*/} >> ./result/assembly_stats_mini.tsv
# done

# for file in ./result/*.fasta
# do
#     echo "======" $file "======\n" >> ./result/assembly_stats.dat
#     assembly-stats $file >> ./result/assembly_stats.dat
#     echo "" >> ./result/assembly_stats.dat
#     assembly-stat $file ${file##*/} >> ./result/assembly_stats_mini.tsv
# done
