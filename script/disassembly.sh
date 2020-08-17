#!/bin/bash
set -ue
REFERENCE=$1
READ=$2
OUTPUT=$3
MIN_CLUSTER=$4
LIMIT=$5
CORES=$6
ROOT=${PWD}
PATH="${PWD}/target/release/:${PATH}"

# ---- Clean up -----
if [ -d ${OUTPUT} ]
then
    rm -r ${OUTPUT}
fi

# ----- Prediction ------
mmmm decompose \
     --output ${OUTPUT} \
     --reads ${READ} --contigs ${REFERENCE} \
     --cluster_num ${MIN_CLUSTER} --threads ${CORES} \
     --limit ${LIMIT}\
     -vv

mkdir -p ${OUTPUT}/no_merge
mmmm decompose \
     --output ${OUTPUT}/no_merge \
     --reads ${READ} --contigs ${REFERENCE} \
     --threads ${CORES} \
     --no_merge

# ---- Assembly(by Flye) ----
set +e
mkdir -p ${OUTPUT}/assemblies/
for reads in $( find ${OUTPUT} -maxdepth 1 -name "*.fasta" )
do
    echo "assembling ${reads}"
    ASM_PATH=${reads%%.fasta}
    INDEX=${ASM_PATH##*/}
    mkdir -p ${OUTPUT}/assemblies/${INDEX}
    genomesize=$(estimate_genome_size ${reads} ${REFERENCE} ${OUTPUT}/assemblies/${INDEX}/temp.fa)
    flye \
        --pacbio-raw \
        ${OUTPUT}/assemblies/${INDEX}/temp.fa \
	    --genome-size ${genomesize} \
        --threads ${CORES} \
        --iterations 10 \
        --out-dir ${OUTPUT}/assemblies/${INDEX}
    rm ${OUTPUT}/assemblies/${INDEX}/temp.fa
done

set -ue
# ---- Align back all reads ----
mkdir -p ${OUTPUT}/last_db
collect_contigs ${OUTPUT}/assemblies/ ${OUTPUT}
for contigs in  $( find ${OUTPUT} -maxdepth 1 -name "*contigs.fasta" )
do
    reads=${contigs%%.contigs.fasta}.fasta
    cd ${OUTPUT}/last_db
    lastdb -R00 temp ${contigs}
    last-train -P${CORES} -Q0 temp ${reads} > temp.matrix
    name=${RANDOM}
    lastal -f maf -P${CORES} -R00 -Q0 -p temp.matrix temp ${reads} \
           > temp${name}.maf
    # Remove unnessary sequence.
    cd ${ROOT}
    remove_low_coverage \
          ${contigs} ${OUTPUT}/last_db/temp${name}.maf > ${OUTPUT}/last_db/temp.fa
    cd ${OUTPUT}/last_db
    mv temp.fa ${contigs}
    lastdb -R00 temp ${contigs}
    last-train -P${CORES} -Q0 temp ${reads} > temp.matrix 
    lastal -f maf -P${CORES} -R00 -Q0 -p temp.matrix temp ${reads} |\
        last-split |\
        maf-convert tab --join 1000 > ${reads%%.fasta}.reads.aln.tab
    last-train -P${CORES} -Q0 reference ${contigs} > temp.matrix
    if [ $? -eq 0 ]
    then
        lastal -f tab -P ${CORES} -R00 -p temp.matrix -Q0 reference ${contigs} \
               > ${contigs%%.fasta}.aln.tab
    else
        lastal -f tab -P ${CORES} -R00 -Q0 reference ${contigs} \
               > ${contigs%%.fasta}.aln.tab
    fi
    cd ${ROOT}
done

cat ${OUTPUT}/*.reads.aln.tab > ${OUTPUT}/allreads.aln.tab
cat ${OUTPUT}/*.contigs.aln.tab > ${OUTPUT}/allcontigs.aln.tab
rm ${OUTPUT}/*.reads.aln.tab ${OUTPUT}/*.contigs.aln.tab


# ---- Create viewer files -----
mmmm create_viewer --assignments ${OUTPUT}/readlist.tsv \
     --contig_aln ${OUTPUT}/allcontigs.aln.tab \
     --contigs ${OUTPUT}/multipartite.fasta \
     --output_dir ${OUTPUT}/viewer/ \
     --read_aln ${OUTPUT}/allreads.aln.tab \
     --reads ${READ} \
     --reference ${REFERENCE}\
     --min_align_length 1000
