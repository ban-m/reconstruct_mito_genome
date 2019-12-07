library("tidyverse")
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}


dataset_diff <- read_tsv("./result/coverage_and_likelihood_diff.tsv")
dataset <- read_tsv("./result/coverage_and_likelihood2.tsv")

g <- dataset_diff %>%
    gather(key = Type, value = Likelihood, -Coverage) %>% 
    filter(Likelihood > -200) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood, color = Type), alpha = 0.3)
generalplot(g,"coverage_and_likelihood_point")

g2 <- dataset %>% filter(Likelihood > -200) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood, color = Type), alpha = 0.3)

g2 <- dataset %>% filter(Coverage > 1) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood), alpha = 0.3) +
    facet_wrap(.~ Type)

g <- dataset_diff %>%
    gather(key = Type, value = Likelihood, -Coverage) %>% 
    filter(Likelihood > -200) %>%
    ggplot() + geom_smooth(aes(x = Coverage, y = Likelihood, color = Type)) 

g2 <- dataset %>% filter(Likelihood > -200) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood, color = Type))


summaries  <- dataset %>% filter(Coverage>1) %>%
    filter(LikelihoodRatio < 100) %>%
    nest(-Seed,-Coverage) %>%
    mutate(data = map(data,function(x) x %>% summarize(mean =mean(LikelihoodRatio)))) %>%
    unnest()

g <- summaries %>% ggplot() + geom_point(aes(x= Coverage,y = mean)) + ## , color = factor(Seed))) +
    stat_function(fun=function(x)exp(-0.22 * x  + 3.58))

g <- dataset %>%
    filter(Coverage>1) %>%
    filter(LikelihoodRatio < 100) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = LikelihoodRatio,color = Seed))

tempdataset  <- dataset %>% filter(Coverage>1) %>%
    filter(LikelihoodRatio < 100)
objective_function <- function(x){
    a <- x[1]
    b <- x[2]
    tempdataset %>%  mutate(residual = (LikelihoodRatio - exp(a*Coverage+b)) ** 2) %>%
        pull(residual) %>% sum()
}
result <- optim(par=c(-0.24,3.6),fn =  objective_function)

g <- summaries %>% ggplot() + geom_point(aes(x= Coverage,y = mean)) + ## , color = factor(Seed))) +
    stat_function(fun=function(x)exp(result$par[1] * x  + result$par[2]))




g3 <-  dataset_diff %>%
    mutate(diff = Correct - Wrong) %>% filter(abs(diff) < 100) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = diff))
generalplot(g,"coverage_and_likelihood_diff")

g <- dataset_diff %>% nest(-Coverage) %>%
    mutate(data = map(data,function(x) x %>% filter(Correct - Wrong < 0) %>% count())) %>%
    unnest() %>%
    ggplot()  + geom_point(aes(x = Coverage, y = n ))




g <- dataset %>% nest(-Coverage,-Type) %>%
    mutate(data = map(data,function(x)summarize(x,mean = mean(Likelihood))))%>%
    unnest() %>%
    ggplot() + geom_point(aes(x = Coverage, y = mean,color = Type))
generalplot(g,"coverage_and_likelihood_mean")

