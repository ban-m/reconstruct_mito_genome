library("tidyverse")
data <- read_tsv("./result/maf.tsv", col_names =FALSE)
g <-data %>%  ggplot() +
    geom_histogram(mapping = aes(x = X2, y = ..density..), bins=100) +
    stat_function(fun = function(x){dnorm(x,mean=mean(data$X2),sd=sd(data$X2))})
ggsave(filename = "./pics/minor_allel_freq.png", g)
g <-data %>%  ggplot() +
    geom_point(mapping = aes(x = X1, y = X2))
ggsave(filename = "./pics/minor_allel_freq_point.png", g)

filtered <- data %>%
    filter((4956 < X1 & X1 < 37653) |
           (109386 < X1 & X1 < 144786) |
           (341256 < X1 & X1 < 367807) )
(37653 - 4956 + 144786 - 109386 + 367807 - 341256)/367807
print(mean(data$X2) + 3*sd(data$X2))
big_maf <- data %>% filter(X2 > mean(data$X2) + 3*sd(data$X2))
big_maf_count <- big_maf %>% count() %>% pull(n)
total <- data %>% count() %>% pull(n)
print(big_maf_count)
print(total)
print(big_maf_count/total * 100)