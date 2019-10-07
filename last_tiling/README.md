# MAF tiling

Author: Bansho Masutani banmasutani@gmail.com

# Desc.

This crate consits of tiling a given reads into a "encoded" version of the read.

# Synopsis

To encode reads into array of units by maf alignment,
```bash
cargo run --release --bin encode -- [read<FASTA>] [maf<MAF>] > [encoded<Binary>]
```

The output file is a binary file encoded via MassegePack encoding.

To see more detail about the program, see documentation.

## ToDo
