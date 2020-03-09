library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
dataset <- read_tsv(args[1], col_names=FALSE)
print(dataset)
g <- dataset %>% rename(Iteration = X2, LogLikelihood =  X3) %>%
    mutate(Replicate = factor(X1)) %>%
    ggplot() + geom_line(aes(x = Iteration, y = LogLikelihood, color = Replicate))
ggsave(paste0("./png/", args[2], ".png"), g)
ggsave(paste0("./pdf/", args[2], ".pdf"), g)
