library("tidyverse")

maxk <- 10
error_rate <- 0.15
dataset <- tibble(k = 2:maxk %>% rep(.,times = . + 1),
       error = unlist(2:maxk %>% sapply(FUN=function(x)seq(0,x))),
       prob = choose(k,error)*(error_rate ** error) * ( (1-error_rate) ** (k-error)))


flow <- function(k,error){
    res <- 0:k %>% (function(i)(1-error)**(k-i) * (error/3)**i)
    names(res) <- 0:k
    res
}
