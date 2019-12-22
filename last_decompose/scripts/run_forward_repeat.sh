ROOT=/grid/ban-m/arabidopsis_thaliana/sequel/assemble/forward_repeat
qsub ./scripts/run.sh ${ROOT}/filtered_read.fasta \
     ${ROOT}/mito.fa \
     ${ROOT}/last_db/initial.tab \
     ${ROOT}/last_db/self.tab \
     forward_repeat

