lib <- .libPaths()
.libPaths(lib[c(1,3,2)])
library("stringi")
library("tidyverse")
source("~/work/generalplot.R")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback/guided_asm_sequel_canu.contigs.fasta.mapback.tsv")

outname <- "genomic_region_coverage"

### ====== Open data =====
data1 <- read_tsv(args[1], col_names=FALSE) %>% mutate(seq = "Sequel")
data2 <- read_tsv(args[2], col_names = FALSE) %>% mutate(seq = "MinION")

data <- bind_rows(data1,data2) %>% rename(chr = X1, position = X2, depth = X3)

g <- data %>% ggplot(mapping = aes(x = position,
                                   y = depth)) + geom_line() +
    facet_grid(seq ~. ) +
    ylim(c(0,250))

generalplot(g,"trim_genomic_region")
