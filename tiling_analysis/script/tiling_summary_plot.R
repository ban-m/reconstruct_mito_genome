library("tidyverse")
source("~/work/generalplot.R")
loadNamespace("cowplot")

### ===== load data ====
args <- commandArgs(trailingOnly =  TRUE)
args <- c("./result/tiling_summary.tsv")

tiling_summary <- read_tsv(args[1])

g <- tiling_summary %>% filter(Type == "GAP") %>% ggplot(mapping = aes(x = Count)) + geom_histogram(bins = 100)  

generalplot(g, "gap_histogram_of_canu")


g <- tiling_summary %>% filter(Type != "GAP") %>% ggplot(mapping = aes(x = Unit, y = Count)) + geom_line() + facet_wrap( . ~ Type, scales = 'free')

name <- "unit_summary_of_canu"
cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                plot = g + cowplot::theme_cowplot(),
                width = 21, height = 7)
cowplot::ggsave(filename = paste0("./png/",name,".png"),
                width = 21,
                height = 7,
                dpi= 320,
                plot = g + cowplot::theme_cowplot())

