# Tiling by kmercount

Author: Bansho Masutani

Email:ban-m@g.ecc.u-tokyo.ac.jp


## Description

Tiling reads with kmers -- in other words, to convert each k-mer into the occurrence count of its "nearest" kmer.

## Synopsis
```bash
cargo run --releae --bin tiling -- [fastq file] [kmer database] > [fastq like file]
```

If this is too huge(it is often the case), one can "smoothing" the kmer counts by given windowsize.
```bash
cargo run --release --bin tiling_and_smoothing -- [fastq file] [kmer database] [window size] > [fastq like file]
```
Note that the [window size] should be greater than the [k].

