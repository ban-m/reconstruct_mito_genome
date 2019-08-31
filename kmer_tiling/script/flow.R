library("tidyverse")
args <- commandArgs(trailingOnly = TRUE)
error <- parse_number(args[1])
k <- parse_integer(args[2])
r <- parse_integer(args[3])

error <- 0.15
k <- 15
r <- 100

flow <- 0:r %>% (function(i){((1-error) ** (k-i)) *( ( error / 3 ) ** i)})
total <- 0:r %>% (function(i){(choose(k,i)) * ((1-error)**(k-i)) * (error ** i)})

data <- tibble(flow = flow, total = total, index = 0:r, cumsum = cumsum(total))
