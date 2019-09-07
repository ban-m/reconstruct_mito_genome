library("tidyverse")
args <- commandArgs(trailingOnly = TRUE)
data <- read_tsv(args[1],col_names = FALSE) %>% rename(count = X2, entropy = X3)

g <- data %>% ggplot(mapping = aes(x = count, y = entropy)) + geom_point()
ggsave(g, args[2])
