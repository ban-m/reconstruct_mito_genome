# Error Correction

Author:Bansho Masutani


# Desc.

To correct errors.

# Synopsis


## ./script/MECAT2.job

A wrapper script for MECAT2.

## diagonal

Exectute the self-alignment between the raw read and corrected version.
```bash
cargo run --release --bin diagonal [RawRead<fastq>] [CorrectedRead<fastq>].fa > ./[Result<mapping>]
```


