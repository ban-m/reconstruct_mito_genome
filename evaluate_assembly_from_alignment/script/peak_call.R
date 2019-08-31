library("tidyverse")
coverages <- read_tsv("./result/mitochondria/coverage.tsv",col_names=FALSE) %>%
    rename(idx = X1, depth =X2)
g <- coverages %>% ggplot(aes(x = idx, y = depth)) + geom_line()

smoothing <- function(xs){
    w <- 200
    len <- length(xs)
    sapply(1:len,function(i) {
        start <- max(i-w,0)
        end <- min(i+w,len)
        mean(xs[start:end])
    })
}

test <- coverages %>% add_column(smoothed=smoothing(smoothing(smoothing(coverages$depth))))
g <- test %>% ggplot(aes(x = idx, y = depth)) + geom_line()

chr1 <- read_tsv("./result/chr1/coverage.tsv",col_names=FALSE) %>%
    rename(idx = X1, depth =X2)
offset <- 100000
g <- chr1 %>% filter(offset < idx & idx < offset + length(coverages$depth)) %>%
    ggplot(aes(x = idx, y = depth)) + geom_line()

combined <- bind_rows(chr1 %>% filter(offset <= idx & idx < offset + length(coverages$depth)) %>%
                      mutate(idx = idx - offset, class = "Chr1"),
                      coverages %>% mutate(class = "Mito"))

g <- combined %>% ggplot(aes(x = idx, y = depth,color = class)) + geom_line() 
