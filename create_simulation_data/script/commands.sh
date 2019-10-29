#!/bin/bash

badread simulate --reference ./data/forward_repeat.fasta --quantity 500x --error_model pacbio --qscore_model pacbio --identity 85,95,3 --length 15000,5000 > ./data/forward_repeat.fastq
badread simulate --reference ./data/reverse_repeat.fasta --quantity 500x --error_model pacbio --qscore_model pacbio --identity 85,95,3 --length 15000,5000 > ./data/reverse_repeat.fastq
