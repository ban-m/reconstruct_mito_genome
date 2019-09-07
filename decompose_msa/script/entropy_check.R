## Entropy file, answer file, output figure path.

library("tidyverse")
args <- commandArgs(trailingOnly = TRUE)
entropy <- read_csv(args[1], col_names = FALSE)
answer <- read_tsv(args[2], col_names= FALSE)


sens <- function(entropy, answer, x){
    tp <- entropy[answer$X1,] %>% filter( x < X2) %>% count() %>% pull(n)
    tp / dim(answer)[1]
}


prec <- function(entropy, answer, x){
    tp <- entropy[answer$X1,] %>% filter( x < X2) %>% count() %>% pull(n)
    tp / entropy %>% filter( x < X2) %>% count() %>% pull(n)
}

t = seq(from=min(entropy$X2),to=max(entropy$X2),length=100)
sensitivity <- sapply(X=t, FUN= function(x){sens(entropy, answer,x)})
precision <- sapply(X=t, FUN=function(x){prec(entropy, answer, x)})

result <- tibble(t = t, precision = precision, sensitivity = sensitivity)
g <- result %>% ggplot(mapping = aes(x = precision, y = sensitivity)) + geom_line()
ggsave(args[3],g)
