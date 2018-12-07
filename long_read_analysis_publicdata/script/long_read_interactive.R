library("tidyverse")

sequel.cov <- read_tsv("/data/ban-m/a_thaliana/mitochondria/long_read/sequel_coverage.wig",col_names = FALSE, col_types = 'cii')
names(sequel.cov) <- c("name","index","cov")

g <- sequel.cov %>% ggplot(aes(x = cov)) + geom_histogram(bins=100)

sequel.cov.mapq10 <- read_tsv("/data/ban-m/a_thaliana/mitochondria/long_read/sequel_coverage_mapq10.wig",col_names = FALSE, col_types = 'cii')
names(sequel.cov.mapq10) <- c("name","index","cov")

g <- sequel.cov.mapq10 %>% ggplot(aes(x=cov)) + geom_histogram(bins=100)

sequel.cov.1 <- read_tsv("/data/ban-m/a_thaliana/mitochondria/long_read/sequel_coverage_1.wig",col_names = FALSE, col_types = 'cii')
names(sequel.cov.1) <- c("name","index","cov")
