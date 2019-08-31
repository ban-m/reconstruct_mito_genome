library("tidyverse")
library("stringi")
source("~/work/generalplot.R")
name <- "16766"
files <- c("11988","16766","1733","18080","2858","3097","3578","6199","6736","7820","7829","8464","9016","9649","9699")
plot_error_rate <- function(name){
    error_rate <- read_tsv("./result/error_rate/" %s+% name,col_names=FALSE) %>%
        rename(index = X1,
               subst = X2,
               ins = X4,
               del = X3,
               depth = X5) %>%
        mutate(depth = depth - del)
    len <- length(error_rate$index)
    formatted <- tibble(index = rep(error_rate$index,4),
                        count = c(error_rate$subst,error_rate$ins,error_rate$del,error_rate$depth),
                        type = c(rep("Subst",len),rep("Ins",len),rep("Del",len),rep("Depth",len)))
    g <-formatted %>%  ggplot(mapping = aes(x = index, y = count, color = type)) + geom_line()
    generalplot(g = g, name = name %s+% "_raw_count")
    g <- formatted %>% ggplot(aes(x = count)) + geom_histogram(bins = 60) + facet_wrap(. ~ type)
    generalplot(g=g, name = name %s+% "_raw_hist")
    error_rate <- read_tsv("./result/error_rate/" %s+% name,col_names=FALSE) %>%
        rename(index = X1,
               subst = X2,
               ins = X3,
               del = X4,
               depth = X5) %>%
        mutate(subst = subst / depth,
               ins = ins/depth,
               del = del/depth)
    len <- length(error_rate$index)

    formatted <- tibble(index = rep(error_rate$index,3),
                        fraction = c(error_rate$subst,error_rate$ins,error_rate$del),
                        type = c(rep("Subst",len),rep("Ins",len),rep("Del",len)))
    g <-formatted %>%  ggplot(mapping = aes(x = index, y = fraction, color = type)) + geom_line()
    generalplot(g = g, name = name %s+% "_freqency")
    g <- formatted %>% filter(fraction < 1) %>% ggplot(aes(x = fraction)) + geom_histogram(bins = 60) + facet_grid(type ~ .)
    generalplot(g=g, name = name %s+% "_raw_hist")
    
}

sapply(files,FUN= plot_error_rate)
