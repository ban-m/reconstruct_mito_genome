
# Synopsis
## upper bound of P-value of the nucleotide conposition at each position
```bash
cargo run --release --bin calc_p_value_at_each_position -- [bam file] [substitution rate] > [tsv file]
```
where [substitution rate] is the rate at which a base would be read as one of other three base by a sequence mistakenly(e.g., around 5% for PacBio, 6% for ONT)

The result is a tab-separated file each line containing tid, position, and its upper bound of p-value.