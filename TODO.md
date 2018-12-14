# Arabidopsis mitochondria

Author:Bansho Masutani
Date:2018/12/07
Email:ban-m@g.ecc.u-tokyo.ac.jp

## Current TODO

- Mitochondria reference guided assembly
  - nanopore
  - sequel

- Read coverage assembly
  - nanopore
  - sequel

- Cycle detection assembly
  - This is what I want to do.
  - No alignment, just ancharing.
  


- Kmer frequency of long reads
  - maybe we can detect...
  - falsy overlap
  - skqerker

- Confirm the badness of current reference
  - nanopore/sequel coverage Mito VS Chr1.
    - pay attention to short reads. We may remove any reads shorter than 1Kb?
    - Also, pay attention to dupliacate mapping. MAPQ should be treated adequetly.
  - Read behavior
    - It seems that, by observing the read orientation and mappingness, it is hard to argue that there are some "child circles" in mito DNA.

## Done

- Read Filtering(PacBio)
  - This *should not* based on alignmnet status.
    - It may causes some bias toward "current" genome, which is highly suspicious.
    - And especially the mitochondrial master circle, it is definetly wrong.
  - By selecting a read among the reads which are from the same well in the PacBio's sequel sequencer.
  - Easy

- de novo assembly and assembly graph visualization
  - Use at least two assembler to construct a assembly graph
    - Racon(minimap2 + rala + racon*2)
    - Canu
    - Wtdbg2
  - Targzzing all assemble into one big file, send to local machine.



- Quality control(Both)
  - Triming suspicious reads
  - It would be O.K. to use existing tools.
  - Easy

- de novo assembly and assembly graph visualization
  - Graph visualizer
    - Bandate in the local computer
