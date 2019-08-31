library("tidyverse")
source("~/work/generalplot.R")
args <- commandArgs(trailingOnly = TRUE)

## -------------

kmerfreq <- read_delim(args[1],delim = " ",col_names = FALSE)
g <- kmerfreq[-nrow(kmerfreq),] %>% ggplot(mapping = aes(x = X1, y = X2)) +
    geom_point() +
    labs(title = args[2])
generalplot(g = g, name = args[2])




