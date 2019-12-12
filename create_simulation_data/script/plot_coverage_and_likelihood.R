library("tidyverse")
dataset <- read_tsv("./result/coverage_and_likelihood.tsv",col_names=FALSE)

tempdata <- dataset %>% filter(X2 > -1000) 

objective_function <- function(x){
    a <- x[1]
    b <- x[2]
    c <- x[3]
    diff <- tempdata %>% mutate(diff = X2 + exp(a*X1 + b) - c) %>% pull(diff)
    sum(diff ** 2)
}

initial_param <- c(-0.24, 3.7, -105)


optim(par = initial_param, fn = objective_function)


## $par
## [1]   -0.2217519    3.8889434 -113.5967033

## $value
## [1] 3707086711

## $counts
## function gradient 
##      142       NA 

## $convergence
## [1] 0

## $message
## NULL
