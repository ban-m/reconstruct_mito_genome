library("tidyverse")
source("~/work/generalplot.R")

args <- commandArgs(trailingOnly = TRUE)

p_values <- read_tsv(args[1]) %>% 
    mutate(minus_log_p = -1*log(p_value)) 

false_positive_rate <- 0.01
corrected_fpr <- 0.01 / p_values %>% count() %>% pull(n)
line <- -1*log(corrected_fpr)

ok <- p_values %>% filter(p_value < corrected_fpr) %>% count() %>% pull(n)
total <-  p_values %>% count() %>% pull(n)

g <- p_values %>%
    mutate(tid = as.character(tid))  %>%  
    ggplot(mapping = aes(x = position , y = minus_log_p, color = tid)) +
    geom_point(size=1) +
    geom_hline(yintercept = line,colour = 'red') +
    labs(title = paste0(ok,"<",corrected_fpr,":out of",total)) +
    facet_grid(tid ~ .)

generalplot(g,args[2])
