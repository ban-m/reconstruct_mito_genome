# Select reads from mitochondria

Author: Bansho Masutani

Email:ban-m@g.ecc.u-tokyo.ac.jp

##

To reconstruct the whole picture of mitochondrial DNA, it is nessesaly to represent each mitochondrial molecule in a given organism(plant). To this end, there are roughly two approach: whole genome assembly and selected assembly.

In the former approach, one assemble, at first, the entire genome of the plant, then determine which contigs are likely to be mitochondrial. After that, all the reads are mapped back to these contigs, filtering out unmaped ones, then the next iteration would be carried out.

This approach, although widely used, have weak points. First, it is still time-consuming to assemble the entire genome. Second, due to the unique structure of mitochondrial genomes in plants, the 'picking up' strategy fails to picking up reads from alternative structure.

The latter approach, namely selected assembly, picking up the reads from mitochondria fristly, then assemble the 'selected' reads. It relies the fact that the mitochondrial DNA is high at copy-number, generating higher coverage than ones from nucleic DNA by tens of several times.

However, some issues still remains. For example, among the genome are transposable elements which are disparsed and numerous.

## Synopsis


### select_read_by_tsv.rs

Pick reads from input fastq. The reads to be extracted is specified by a simple line-by-line ID file.
```
cargo run --release --bin select_read_by_tsv -- [fastq file] [line-by-line ID file]
```
### kmer_freq_vector

Convert a readset into vector of k-mer frequency.
```
cargo run --release --bin kmer_freq_vector -- [fasta file] [K]
```

The output format is as follows:

|ID| AAAA | AAAC | AAAG | ... |
|--|------|------|------|-----|
|m_214324_323_34| 0.0232|0.2321|0.003232|.....|
