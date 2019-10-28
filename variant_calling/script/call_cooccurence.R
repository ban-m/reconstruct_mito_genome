library("tidyverse")
variant_data <- read_tsv("./result/variant_called_sites.tsv")

variant_data_filtered <- variant_data %>%
    filter( share > 2)  %>%
    mutate(mlpvalue = -1 * pbinom(q=share, size = mac2, prob = mac1/tot1, lower.tail=FALSE,log.p=TRUE)) 

g <- variant_data_filtered %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = pvalue))

ggsave(filename="./pics/variant_call_pvalue.png",plot = g)

g <- variant_data %>%
    ggplot() +
    geom_histogram(mapping = aes( x = share), bins = 60)
