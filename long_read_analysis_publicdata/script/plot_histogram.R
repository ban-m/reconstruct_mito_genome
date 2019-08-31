library("tidyverse")
args <- commandArgs(trailingOnly = TRUE)

histogram <- read_delim(file= args[1], delim = " ", col_names = FALSE)
g <- histogram[-nrow(histogram),] %>% ggplot(mapping = aes(x = X1, y = X2)) + geom_point();
ggsave(filename = paste0(args[1],".png"),plot=g)
