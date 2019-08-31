## Synopsis

```
cargo run --release --bin kmer_count -- [path to fastq file] [k] [lower limit] [upper limit] > [output file]
```
output: tsv, [# of count]\t[# of kmers with that count]