# Repeat Analysis

Author: Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>

# Contents

## Summary


## Synopsis

```bash
cargo run --release -- [Reference<Fasta>] [SelfAlignment<LastTAB>] [Locations<TSV>] > [TSV]
```
,where Locations should be like "ID\tLocation A Start\tLocation A End\tLocation B Start\tLocation B End\t".
The output format is

Column | Desc.
-------|-------
  1    |  ID
  2    | Start position of the nearest repeat of A
  3    | End position of the nearest repeat of A
  4    | Start position of the nearest repeat of B
  5    | End position of the nearest repeat of B
  6    | Length of the nearest repeat
  7    | Distance from A to the nearest repeat
  8    | Distance from B to the nearest repeat
-------|-------

Last should be used to annotate the reference genome.

