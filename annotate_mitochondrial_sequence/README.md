# Annontate mitochondrial sequence


This is a rough script and binary to generate a annotation file of given intput,
and plotting a figure.


# Installation
```bash
git clone [this repo]
cargo build --release
```

# Usage
```bash
./script/annotate.sh ${INPUT FASTA} ${OUTPUT_DIR} ${PREFIX} ${DATABASE} ${ANNOTATION_FILE}
```

Usually, the ${DATABASE} and ${ANNOTATION_FILE} would be the ones from NCBI mitochondrial genomes.

# Dependencies

Currently, this scripts uses `tRNAscan-SE` version 2 and `last`.
Please install these programs and put them where the program could find(i.e., ${PATH})