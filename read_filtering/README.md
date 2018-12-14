# Read filtering

## Summary
This crate is to filtering out raw quality or suspicious reads.
This is probably needed to trim false positive, and/or artificial reads produced by ONT nanopore sequencer and pacbio sequel sequencer.

Note that, as a rule of thum, I only fitering out reads which I can conclude useless confidently. It is O.K to leave some false positive reads as is. This is because I can detect these reads in the following workflows.

The procedure employed varies depending on the platform. For more details, see below.


### PacBio

I use only one filter: the uniquness of well. Specifically, in the header of fastq entry, a reads begins with
lines like `> @m54113_160913_184949/4194800/0_1273`. Semantically, it means `@m[experiment number]_[sequencing date]_[sequencing time]/[well number]/[movie start time]_[movie end time]`. So, if one want to expect the coverage is to follow Poisson distribution, one should pick a reads from the reads sharing the same well numbers.

I pick the longest reads from these reads.


### Nanopore

It is probably more difficult than that of sequel. I postpone it.


## Synopsis
```
cargo run --release -- "nanopore"|"pacbio" fastq.fq|fast.fa > [output].(fa|fq)
```
