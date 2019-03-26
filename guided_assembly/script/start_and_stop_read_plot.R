lib <- .libPaths()
.libPaths(lib[c(1,3,2)])
library("stringi")

library("tidyverse")
source("~/work/generalplot.R")
args <- commandArgs(trailingOnly = TRUE)
args <- c("/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback/guided_asm_sequel_canu.contigs.fasta.mapback.tsv",
          "/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback/guided_asm_sequel_canu.contigs.fasta.mapback.peaks.tsv"
)

outname <- basename(tools::file_path_sans_ext(args[1]))

### ====== Open data =====
data <- read_tsv(args[1], col_names=TRUE)
peaks <- read_tsv(args[2],col_names=TRUE) %>%
    mutate(type = ifelse(number_of_start_read < number_of_stop_read, "stop","start"))

g <- data %>% ggplot(mapping = aes(x = position,
                                   y = number_of_start_read)) + geom_point() +
    facet_grid(refname ~. ) +
    ylim(c(0,250))

sorted_data <- bind_rows(data %>% select(-number_of_start_read) %>%
                         rename(number_of_read = number_of_stop_read) %>%
                         mutate(type = "stop"),
                         data %>% select(-number_of_stop_read) %>%
                         rename(number_of_read = number_of_start_read) %>%
                         mutate(type = "start"))

g <- sorted_data %>% ggplot(mapping = aes(y = number_of_read, x = position,
                                          color = type)) +
    geom_line(alpha=0.7) +
    facet_grid(refname ~ . ) + 
    ylim(c(0,250))
generalplot(g,paste0(outname,"start_and_stop"))

g <- sorted_data %>% ggplot(mapping = aes(y = number_of_read, x = position,
                                          color = type)) +
    geom_line(alpha=0.7) +
    geom_segment(mapping = aes(x = position,y = 0, xend = position, yend = 200,
                               alpha = 0.3,
                               size = 10,
                               colour = type),
                 data = peaks) + 
    facet_grid(refname ~ . ) + 
    ylim(c(0,250))
generalplot(g,paste0(outname,"start_and_stop_withpeak"))

#### ==== For Canu ====
## plot_by_name <- function(data, ctg_name){
##     g <- data %>% filter(refname == ctg_name) %>%
##         ggplot(mapping = aes(x = position, y = number_of_read,
##                              color = type)) +
##         geom_line() 
##     generalplot(g,outname %s+% "." %s+% ctg_name %s+% ".start_and_stop")
##     1
## }

## for_canu <- function(df){
##     data %>% pull(refname) %>% unique() %>% 
##         sapply(function(x) plot_by_name(df,x))
## }

## for_canu(sorted_data)
