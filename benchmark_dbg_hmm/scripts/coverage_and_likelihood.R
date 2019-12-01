library("tidyverse")


dataset_diff <- read_tsv("./result/coverage_and_likelihood_diff.tsv")
dataset <- read_tsv("./result/coverage_and_likelihood.tsv")

g <- dataset_diff %>%
    gather(key = Type, value = Likelihood, -Coverage) %>% 
    filter(Likelihood > -200) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood, color = Type), alpha = 0.3) 

g2 <- dataset %>% filter(Likelihood > -200) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood, color = Type), alpha = 0.3)


g <- dataset_diff %>%
    gather(key = Type, value = Likelihood, -Coverage) %>% 
    filter(Likelihood > -200) %>%
    ggplot() + geom_smooth(aes(x = Coverage, y = Likelihood, color = Type)) 

g2 <- dataset %>% filter(Likelihood > -200) %>% 
    ggplot() + geom_smooth(aes(x = Coverage, y = Likelihood, color = Type))


g <-  dataset_diff %>%
    mutate(diff = Correct - Wrong) %>% filter(abs(diff) < 100) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = diff))

dataset_diff %>% nest(-Coverage) %>%
    mutate(data = map(data,function(x) x %>% filter(Correct - Wrong < 0) %>% count())) %>%
    unnest() %>%
    ggplot()  + geom_point(aes(x = Coverage, y = n ))


g <- dataset %>% nest(-Coverage,-Type) %>%
    mutate(data = map(data,function(x)summarize(x,mean = mean(Likelihood))))%>%
    unnest() %>%
    ggplot() + geom_point(aes(x = Coverage, y = mean,color = Type))

g <- 
