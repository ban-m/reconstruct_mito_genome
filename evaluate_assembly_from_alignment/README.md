# Evaluation of mitochondrial assembly

Author: Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>

# Synopsis

```bash
cargo run --release -- [sam file] [query file] [reference file] [output directly]  [contig name]
```
If there are no contig with name [contig name], the list of contig names would appair.

This program would create the output directly, if do not exists, and output some metrics into that directly.

The metrics are as follows:

1. The coverage: # of positions with coverage > 1 / # of positions with coverage = 0
   This metrics serves rather detectes sequencing completeness. In other words, we should see a coverage zero with probability exp(- average depth). Thus, if the coverage is quite small exept large depth, it seems to be something 'contaminate' in the assembly.

2. The monotonicity of coverage: the mean of the coverage and the variance.
   This metrics serves sequencing monotonicity. Roughly speaking, the coverage would follow a Poisson distribution with parameter \lambda. Thus, the mean should be equal to the variance, with some disturbance. If it difference by far, something wrong. For example, there should be duplications, deletions, and/or recommbinations.

3. Read clipping: Total mapped reads, # of head/tail clip, # of doulble clip, and # of junction clips.
4. Error profiles at each position. In other words, # of mutation, # of insertion after that location and # of deletion at that position.