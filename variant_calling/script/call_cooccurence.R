library("tidyverse")
variant_data <- read_tsv("./result/variant_called_sites.tsv")

variant_data_filtered <- variant_data %>%
    filter( share > 2)  %>%
    mutate(mlpvalue = -phyper(q=share, m=mac1,n=tot1-mac1, k = mac2, lower.tail=FALSE, log.p=TRUE))  %>%
    filter(!is.na(mlpvalue))


temp <- variant_data_filtered %>% filter(mlpvalue > 10) %>% select(pos1,pos2)
unique_sites <- c(temp$pos1, temp$pos2) %>%
    sort() %>% unique() %>% length()
rm(temp)

g <- variant_data_filtered %>% filter(mlpvalue < 30) %>% 
    ggplot() +
    geom_histogram(mapping = aes(x = mlpvalue), bins = 100)
ggsave(filename="./pics/variant_call_pvalue.png",plot = g)


g <- variant_data_filtered %>%
    filter(10 < mlpvalue & mlpvalue < 30) %>%
    mutate(dist = abs(pos1-pos2)) %>%
    select(dist) %>%
    ggplot() +
    geom_histogram(mapping = aes(x = dist))
ggsave(filename="./pics/variant_call_dist.png",plot=g)

g <- variant_data_filtered %>% filter(10 < mlpvalue & mlpvalue < 30) %>%
    sample_frac(0.01) %>% 
    mutate(size = mlpvalue / 10) %>% 
    ggplot() +
    geom_point(mapping = aes(x=pos1,y = pos2, size = size))
ggsave(filename = "./pics/variant_call_pvalue_point.png", plot = g)

g <- variant_data_filtered %>%
    filter(mlpvalue == Inf) %>% 
    ggplot() +
    geom_point(mapping = aes(x=pos1,y = pos2))
ggsave(filename = "./pics/variant_call_pvalue_point_inf.png", plot = g)

variant_data_filtered %>%
    select(pos1,pos2,mlpvalue) %>%
    filter(mlpvalue > 9) %>%
    mutate(mlpvalue = ifelse(mlpvalue == Inf, 10000, mlpvalue)) %>% 
    write_tsv("./result/variant_cooccurence.tsv")
