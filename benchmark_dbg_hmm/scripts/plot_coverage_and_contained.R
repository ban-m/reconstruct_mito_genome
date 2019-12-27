library("tidyverse")
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}
dataset <- read_tsv("./result/coverage_and_contained_self.tsv")
g <- dataset %>% filter(Coverage > 1) %>% sample_frac(0.1) %>%
    gather(key = Type, value = LK, -Coverage, -Seed) %>%
    ggplot() + geom_point(aes(x=Coverage, y = LK, color= Type),alpha=0.3)
generalplot(g,"coverage_and_contained_overview")


g <- dataset %>% filter(Coverage > 1) %>% sample_frac(0.1) %>%
    gather(key = Type, value = LK, -Coverage, -Seed) %>%
    ggplot() + geom_point(aes(x=Coverage, y = LK)) + facet_wrap(.~Type)


tempdf <- dataset %>%
    mutate(diff = ContainedSelf - NotContained) %>%
    filter(diff > 0) %>% mutate(diff = log(diff)) %>% select(-Seed)


linear_reg <- lm(diff~Coverage, data= tempdf)

objective_function <- function(x){
    a <- x[1]
    b <- x[2]
    dataset %>%  mutate(diff = ContainedSelf - NotContained) %>%
        mutate(residual = (diff - exp(a*Coverage+b)) ** 2) %>%
        pull(residual) %>% sum()
}

init_param <- c(linear_reg$coefficients[2], linear_reg$coefficients[1])
result <- optim(par=c(-0.24,3.6),fn =  objective_function)

param <- result$par
g <- dataset %>% filter(Coverage > 1) %>% sample_frac(0.1) %>%
    mutate(ContainedSelf = ContainedSelf - exp(param[1]*Coverage + param[2])) %>% 
    gather(key = Type, value = LK, -Coverage, -Seed) %>%
    ggplot() + geom_point(aes(x=Coverage, y = LK, color= Type),alpha=0.3)
generalplot(g, "coverage_and_contained_calib")


g <- dataset %>% filter(Coverage > 1) %>% sample_frac(0.1) %>%
    mutate(ContainedSelf = ContainedSelf - exp(param[1]*Coverage + param[2])) %>% 
    gather(key = Type, value = LK, -Coverage, -Seed) %>%
    ggplot() + geom_point(aes(x=Coverage, y = LK)) + facet_wrap(.~Type)
