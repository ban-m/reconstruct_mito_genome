
#!/bin/bash
set -ue
ROOT=${PWD}
READS=${PWD}/../create_simulation_data/data/complex/read_complex1.fa
OUTPATH=${PWD}/result/complex1
REFERENCE=${PWD}/../create_simulation_data/data/complex/reference.fa
ALIGNMENT=${PWD}/result/complex1/last_db/reads.tab
SELF=${PWD}/result/complex1/last_db/self.tab
ASSIGN=${PWD}/result/complex1/yg/readlist3.tsv
cargo run --release --bin create_viewer --\
      ${READS} ${REFERENCE} ${ALIGNMENT} ${SELF} ${OUTPATH} ${ASSIGN}

