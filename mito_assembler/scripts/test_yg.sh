#!/bin/bash
set -ue
ROOT=${PWD}
READS=${PWD}/../create_simulation_data/data/mock_genome_read.fa
OUTPATH=${PWD}/result/mock_genome/
REFERENCE=${PWD}/../create_simulation_data/data/mock_genome_ref.fa
qsub ./scripts/yellow_gate.job ${READS} ${OUTPATH} ${REFERENCE}
