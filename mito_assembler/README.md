# Mito Aassembler

Author: Bansho Masutani
Email: ban-m@g.ecc.u-tokyo.ac.jp

# Desc.

This is a tiny "disassembler" of plant mitochondial genome.

# Requirements

- Rust lang
- SMRTPipe version 6.0.
- Fastq file for target organism

- Docker or the softwares below:
  - minimap2
  - LAST
  - Flye


# Directry structures


# Contents

## Program description

### Rust

- filter_genomic_reads
Input: Read<Fasta>
Output: Read<Fasta>
Requirements: Alignments<Paf> from stdin.
Filtering out reads mapped to genomic regions.

- select_mito_reads.rs
Input: Read<Fasta>
Output: Read<Fasta>
Requirements: Alignments<PAF> from stdin
Select all the reads mapped to mitochondrial genome.

- split_subreads.rs
Input: None
Output: Read<Fasta>
Requirements: Read<Fasta> should be supplied from stdin. The ID of the Input should follow the PacBio's nomenclature.
Select one of a read from the same well.

- filter_low_quality.rs
Input: None
Output: Read<Fasta>
Requirements: Read<Fasta> should be supplied from stdin.
Filter out reads with 6-mer entropy less than 9.

- clip_self_chimera.rs
Input: None
Output: Read<Fasta>
Requirements: Read<Fasta> should be supplied from stdin.
Chop input reads at adaptor chimera.

- split_repeat.rs
Input: Contigs<Fasta>, and self-vs-self alignment<LastTAB>
Output: Fasta
Collupse all the repeat and output the remmaining contigs, with the split repeats.

- entropy.rs
Input: Read<Fasta>
Output: Entropies<TSV>
Requirements: No
Calcurate the k-mer entorpy for each read for k = 2,3,4,5,6. Read information would be lost.

- decompose.rs
Input: Read<Fasta>, Alignment<LASTTAB>, Contigs<Fasta> , self-self alignment<LastTAB>, and the location where the output would be written<Path>,
Output: None,
Requirements: None
Decompose the reads into blocks.

- calc_error_profile.rs
Input: Alignemnts<LASTTAB>, Read<Fasta>, and Contigs<Fasta>
Output: ErrorProfile<TXT>
Calculate error rates of given dataset. The alignments should be alignments from the reads the contigs.

### ./scripts

- workflow.job
Main script. Invoked by run.sh.

- run.sh
Wrapper for workflow.job.