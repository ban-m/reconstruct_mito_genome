library("tidyverse")
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}


dataset_diff <- read_tsv("./result/coverage_and_likelihood_diff.tsv")
dataset <- read_tsv("./result/coverage_and_likelihood.tsv")

g2 <- dataset %>% filter(Coverage > 1) %>% filter(LikelihoodRatio < 100) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = LikelihoodRatio), alpha = 0.3)
generalplot(g2,"coverage_and_likelihood_ratio_point")
g2 <- dataset %>% filter(Coverage > 1) %>% filter(LikelihoodRatio < 100) %>% 
    ggplot() +
    geom_smooth(aes(x = Coverage, y = LikelihoodRatio, color = factor(Seed)), alpha = 0.3)
generalplot(g2,"coverage_and_likelihood_ratio_smooth")


summaries  <- dataset %>% filter(Coverage>1) %>%
    filter(LikelihoodRatio < 100) %>%
    nest(-Seed,-Coverage) %>%
    mutate(data = map(data,function(x) x %>% summarize(mean =mean(LikelihoodRatio)))) %>%
    unnest()

tempdataset  <- dataset %>% filter(Coverage>2) %>%
    filter(OrigignalLK > -200) 

linear_reg <- lm(log(LikelihoodRatio) ~ Coverage,data=tempdataset)
objective_function <- function(x){
    a <- x[1]
    b <- x[2]
    tempdataset %>%  mutate(residual = (LikelihoodRatio - exp(a*Coverage+b)) ** 2) %>%
        pull(residual) %>% sum()
}
init_param <- c(linear_reg$coefficients[2], linear_reg$coefficients[1])
result <- optim(par=c(-0.24,3.6),fn =  objective_function)

g <- dataset %>% filter(Coverage > 1 ) %>% filter(LikelihoodRatio < 100) %>% 
    ggplot() + geom_point(aes(x= Coverage,y = LikelihoodRatio, color = factor(Seed))) +
    stat_function(fun=function(x)exp(result$par[1] * x  + result$par[2]))
generalplot(g, "coverage_and_likelihood_ratio_with_regress")

a <- result$par[1]
b <- result$par[2]
g <- dataset %>% filter(Coverage > 1 ) %>% filter(LikelihoodRatio < 100) %>%
    mutate(LikelihoodRatio = LikelihoodRatio - exp(a*Coverage+b)) %>% 
    ggplot() + geom_point(aes(x= Coverage,y = LikelihoodRatio, color = factor(Seed)))


generalplot(g, "coverage_and_likelihood_ratio_sub_regress")
g <- dataset %>% filter(Coverage > 1 ) %>% filter(LikelihoodRatio < 100) %>%
    mutate(LikelihoodRatio = LikelihoodRatio - exp(a*Coverage+b)) %>% 
    ggplot() + geom_smooth(aes(x= Coverage,y = LikelihoodRatio, color = factor(Seed)))
generalplot(g, "coverage_and_likelihood_ratio_sub_regress_smooth")



g <- dataset_diff %>%
    gather(key = Type, value = Likelihood, -Coverage) %>% 
    filter(Likelihood > -200) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Likelihood, color = Type), alpha = 0.3)
generalplot(g,"coverage_and_likelihood_point")

g <- dataset_diff %>%
    gather(key = Type, value = Likelihood, -Coverage) %>% 
    filter(Likelihood > -200) %>%
    ggplot() + geom_smooth(aes(x = Coverage, y = Likelihood, color = Type))
generalplot(g,"coverage_and_likelihood_smooth")



g3 <-  dataset_diff %>%
    mutate(diff = Correct - Wrong) %>% filter(abs(diff) < 100) %>% 
    ggplot() + geom_point(aes(x = Coverage, y = diff))
generalplot(g,"coverage_and_likelihood_diff")


g <- dataset_diff %>% nest(-Coverage) %>%
    mutate(data = map(data,function(x) x %>% filter(Correct - Wrong < 0) %>% count())) %>%
    unnest() %>%
    ggplot()  + geom_point(aes(x = Coverage, y = n ))


