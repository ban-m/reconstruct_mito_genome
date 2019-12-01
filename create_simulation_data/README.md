# Create simulated dataset


## Contents


### Synopsis


- main.rs
Input: Arabidopsis mitochondria, enhanced genome<FASTA>
Output: None
Requirement: None
Create reference genomes for simulations. It would create `./data/forward_repeat.fasta` and `./data/reverse_repeat.fasta`, which have two or three fasta file correspondig to the recommbination procedure, respectively.

- compare_accreagte_strategy.rs
Input: None
Output: TSV; A record consists of performance of naive aggregation, weighted sum, weighted log likelihood, distance, coverage and the length of the pattern.

- simulated_data.rs
Input: None
Output: TSV; A record consists of performance of weighted sum, distance, coverage, and the length of the pattern.

- create_mock_genome.rs
Input: The length of the genome to be created.
Output: A Fasta record.
Create reference genomes for simulations. It would first create a random genome with specified length, then introduce some mutations. Both depths of these contigs are set to 1.0. These genomes are circular. The error profile is as follows: substitution 0.2%, insertion 0.2%, and deletion 0.2%.


- create_complext_structure.rs
Input: The length of the genome to be created<Int> and the path to the directory to output<Path>
Output:None
Create reference genomes for simulations. It would first create a circular genome with four building blocks, A, B, C, and D. This program writes this genome as "reference.fa". This genome is circular. Depth is 1.0.
Then, this program create another structure, A->D->C. This sequence has linear structure and has depth=0.5. Then, the reference and this structure is output as "complex1.fa".
Also, two structure A->D(Circular, depth=0.5) and C->D(Circular, depth=0.5) would be created and output as "complex2.fa". Each structure has error rate of 0.6% from the reference sequence independently each other.

- predict_mockdata.rs
Input:
Output: TSV.
This binary should be called from ./scripts/mock_genome_workfow.job.

- compare_aggregate_strategy_on_mockdata.rs
Input: Reads<Fasta>, Reference<Fasta>, and alignemnts<LastTAB>
Output: TSV
Make predictions.

- em_algorithm_check.rs
Input:None
Output: TSV
Assess EM algorthm.