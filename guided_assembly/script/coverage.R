lib <- .libPaths()
.libPaths(lib[c(1,3,2)])

library("tidyverse")
library("stringi")
source("~/work/generalplot.R")


### ======= Load data =======
# args <- c("/grid/ban-m/arabidopsis_thaliana/sequel/guided_asm/mapback_after_split/guided_asm_sequel_canu_split.contigs.mapback.wig")
args <- commandArgs(trailingOnly = TRUE)
data <- read_tsv(args[1],col_names = FALSE) %>% rename(contig_name = X1,
                                                       position = X2,
                                                       depth = X3)
print(args[1])
print(summary(data$depth))
filestem <- basename(args[1])

res <- data %>% nest(-contig_name) %>%
    mutate(data = map(data,function(df) df %>% summarize(mean = mean(depth),
                                                         sd = sd(depth)))) %>%
    unnest() %>%
    mutate(name = filestem) %>%
    select(name, contig_name, mean, sd)

write_csv(res, "./data/" %s+% filestem %s+% "_coverage.csv",col_names=FALSE)
### ==== Plot data =======
trim_upper_lower_threshold <- function(df) {
    # Trim upper/lower 3 % of the data.
    u_threshold <- df %>% pull(depth) %>% quantile(0.97)
    l_threshold <- 0
    df %>% filter(depth < u_threshold & l_threshold < depth)
}

trimed_data <- data %>% nest(-contig_name) %>%
    mutate(data = map(data,trim_upper_lower_threshold)) %>% unnest()

longest_contig_name <- trimed_data %>% nest(-contig_name) %>%
    mutate(data = map(data, function(df) df %>% summarize(length = length(depth)))) %>%
    unnest() %>%
    arrange(desc(length)) %>%
    head(n=1) %>%
    pull(contig_name)

g <- trimed_data %>%
    filter(contig_name == longest_contig_name) %>% 
    ggplot(aes(x = position,y = depth)) + geom_line()  + ylim(c(0,NA))

generalplot(g,filestem %s+% ".positionwize")

g <- trimed_data %>%
    ggplot(aes(x = position,y = depth)) + geom_line()  + ylim(c(0,NA)) + facet_grid(contig_name ~ .)
generalplot(g,filestem %s+% ".positionwize.contigwize")

## g <- data %>% ggplot(mapping = aes(x = depth)) + geom_histogram(bins=60) + facet_grid(contig_name ~ . )
## generalplot(g,filestem %s+% ".histogram")

## prev_depth <- c(0,data$depth)[1:dim(data)[1]]
## diff <- data %>% mutate(diff = depth - prev_depth)

## g <- diff %>% ggplot(mapping = aes(x = position,y = diff)) + geom_point() +
##     ylim(c(-1000,1000))
## generalplot(g,filestem %s+% ".diff.positionwize")
## g <- diff %>% filter(abs(diff) < 1000) %>%
##     ggplot(mapping = aes(x = diff)) + geom_histogram(bins=60) 
## generalplot(g,filestem %s+% ".diff.histogram")

#### ===== Only for Canu ====
## plot_by_name <- function(data, ctg_name){
##     g <- data %>% filter(contig_name == ctg_name) %>%
##         ggplot(mapping = aes(x = position, y = depth)) + geom_line()
##     generalplot(g,filestem %s+% "." %s+% ctg_name %s+% ".histogram")
##     1
## }

## for_canu <- function(df){
##     data %>% pull(contig_name) %>% unique() %>% 
##         sapply(function(x) plot_by_name(df,x))
## }

## for_canu(data)


## threshold <- data %>% pull(depth) %>% quantile(.98)
## filtered <- data %>% filter(depth < threshold) %>% arrange(position)
## binwidth <- 1000
## binning <- function(df){
##     length <- df %>% tail(n=1) %>% pull(position)
##     bins <- tibble(position = 1:length,
##                    binnum = position %/% binwidth)
##     right_join(df, bins) %>% nest(-binnum) %>%
##         mutate(data = map(data,function(df) df %>% summarize(mean = mean(depth,na.rm = TRUE),
##                                                              size = sum(!is.na(depth))))) %>%
##         unnest()
## }

## binned <- data %>% nest(-contig_name) %>% mutate(data = map(data,binning)) %>% unnest()
