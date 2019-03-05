l <- .libPaths()
.libPaths(l[c(1,3)])

library("tidyverse")

coverage <- read_tsv("./result/sequel_coverages.tsv",col_names = FALSE) %>%
    rename(id=X1, coverage=X2)
assign <- read_delim(delim=' ',"./result/sequel_mapped_ctg.tsv", col_names = FALSE) %>%
    rename(id=X1, contig=X2)

data <- right_join(coverage, assign)

data <- read_tsv("./result/data.tsv")
g <- data %>% ggplot(mapping=aes(x = contig, y = coverage)) + geom_jitter()

pdf("./pdf/density.pdf")
plot(density(data$1))
dev.off()
