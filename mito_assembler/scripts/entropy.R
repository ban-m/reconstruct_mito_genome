library("tidyverse")
dataset <- read_tsv("./data/entropy.tsv")
g <- dataset %>% ggplot(mapping = aes(x = entropy)) + geom_histogram() + facet_wrap(k ~ . )
plot(g)
