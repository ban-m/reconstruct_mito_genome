library("tidyverse")

p_values <- read_tsv("./result/sequel_canu_P_value.tsv") %>%
    mutate(minus_log_p = -1*log(p_value)) 


false_positive_rate <- 0.01
corrected_fpr <- 0.01 / p_values %>% count() %>% pull(n)
line <- -1*log(corrected_fpr)


g <- p_values %>% 
    ggplot(mapping = aes(x = position , y = minus_log_p)) +
    geom_point () +
    geom_hline(yintercept = line,colour = 'red')
ggsave(filename= "./pdf/sequel_canu_p_values.pdf",g)
