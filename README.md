# Landscape of plant mitochondrial genome via dis-assembly approach.



## Short Description



This repository is for my research on plant mitochondrial genomes, which I call "mitogenome" in this README.md. Shortly and roughly speaking, mitogenomes of plants are different from those of animals or fungus. For example, 

- The sizes are much larger.
- Genomes contain some repetitive regions, and
-  (repeat-)mediated recombinations create divergent structures called "multipartite structure"s.

I'm trying to reconstruct these multipartite structures in the plant mitogenomes. In contrast to usual assembly workflows, I'd call my approach as "dis-"assembly, as I leverage existing reference genomes and "find" other structures by clustering.



## Properties

- Author: Bansho Masutani at University of Tokyo
- Last Update: 2020/05/29
- Contact: ban-m@g.ecc.u-tokyo.ac.jp



## Data Availability

The following list is the list of accessions/URLs to the data used in this research:

- Seven accessions of Arabidopsis thaliana via [Sequence Read Archive]( https://www.ebi.ac.uk/ena/data/view/PRJEB31147)
- One Ler accession from [PacBio's official repository](https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/)
- One Col-0 accession sequenced newly in the National Institute of Genetics in Japan (BioSample: SAMD00224135-SAMD00224139).

All of these data were generated by Sequel Systems. For more details, e.g., the used chemistory, please visit the links.

The Gen Bank IDs of the reference sequences used in this research are BK010421.1, JF729100, and JF729102. BK010421.1 is the reference sequence I'm using generally.

Also, the supplementary plots are available at the following locations. "Circos," "linear," and "no-marge" represent a circos plot (like in the paper), a dot plot between assembled contig and the reference, and a circos plots without any "aggressive merging," respectively.

|Strain Name | Circos | Linear | no_merge |
|:----------:|:------:|:------:|:--------:|
|pacbio | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/pacbio/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/pacbio/no_merge.html)|
|ler | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/ler/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/ler/no_merge.html)|
|col0| [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/col0_1106_exp2/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/col0_1106_exp2/no_merge.html)|
|sha | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/sha/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/sha/no_merge.html)|
|cvi-0 | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/cvi/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/cvi/no_merge.html)|
|an1 | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/an1/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/an1/no_merge.html)|
|c24 | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/c24/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/c24/no_merge.html)|
|kyoto | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/kyo/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/kyo/no_merge.html)|
|eri-1 | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/eri/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/eri/no_merge.html)|
|pacbio_ler | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/pacbio_ler/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/pacbio_ler/no_merge.html)|
|ler_ler | [View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/ler_ler/circos.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly//linear.html)|[View](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/ler_ler/no_merge.html)|

Note: We used BK010421.1 as the reference genome except pacbio_ler and ler_ler where we used JF729100 as the reference genome.

For those interestead in the raw data, I attach a link to the tar.gzed file containing all of the result;

- The result of the synthetic data is [Here](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/simulated_data_result.tar.gz)(1.9G).
- The result of the real data is [Here](https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/mitochondria_assembly/disassembly.tar.gz)(16G).

## Reproducibility


For the synthetic dataset, as BadRead program does not have any options to fix the seed for a quasi-random number generator(2019/05/29), I gnu-zipped the entire dataset as a single tar-ball. [Here is the link]().


To reproduce our result for Fig2 a,b,c, go to ./simulated_data directory and run `./script/local_re_run.sh`. It would compile & run the experiment and generate plots in `png` or `pdf` directory. Note that it would take a few days to complete, as it runs the same program with different seeds again and again.

To reproduce our result for the rest of Fig2
 and Fig3, go disassembly directory and run `./scirpt/local_run.sh`. It might take a week to complete in a local computer, so I highly recommend executing it in parallel. To do that, if your environment's job-schedular is SGE or SGE-compatible, run `./script/real_dataset.sh` and `./script/complex.sh`. If you are using other schedulers, please contact me.



## Reference

- Jiao, W., Schneeberger, K. Chromosome-level assemblies of multiple *Arabidopsis* genomes reveal hotspots of rearrangements with altered evolutionary dynamics.
  *Nat Commun* **11,** 989 (2020). https://doi.org/10.1038/s41467-020-14779-y

