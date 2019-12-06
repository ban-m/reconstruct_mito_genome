#!/bin/bash
set -ue
qsub -sync yes ./script/run_workflow.job \
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/mock_genome_read.fa \
     ${PWD}/result/simulated_genome/ \
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/mock_genome_ref.fa

qsub -sync yes ./script/run_workflow.job\
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/flip_repeats/read_forward_repeat.fa\
     ${PWD}/result/forward_repeat/ \
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/flip_repeats/reference.fasta

qsub -sync yes ./script/run_workflow.job\
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/flip_repeats/read_reverse_repeat.fa\
     ${PWD}/result/reverse_repeat/ \
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/flip_repeats/reference.fasta

qsub -sync yes ./script/run_workflow.job\
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/complex/read_complex1.fa \
     ${PWD}/result/complex1/ \
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/complex/reference.fa

qsub -sync yes ./script/run_workflow.job\
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/complex/read_complex2.fa \
     ${PWD}/result/complex2/ \
     ~/work/arabidopsis_mitochondria/create_simulation_data/data/complex/reference.fa
