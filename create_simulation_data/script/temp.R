

objective_function <- function(x){
    a <- x[1]
    b <- x[2]
    c <- x[3]
    dataset %>% filter(3 < X1) %>%
        mutate(diff = c - (X3 +exp(a*X1 +b))) %>% mutate(diff = abs(diff)) %>%
        pull(diff) %>% sum()
}
result <- optim(par = c(-0.1,3,-100), fn = objective_function)


gf <- dataset %>% filter(X1 > 2) %>% mutate(LK = X3 + exp(1.15-X1 * 0.02)) %>%
    ggplot() + geom_point(aes(x = X1, y = LK))

gp <- dataset %>% filter(X1 > 2) %>% 
    ggplot() + geom_point(aes(x = X1, y = X3))
