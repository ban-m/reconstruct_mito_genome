library("tidyverse")
source("~/work/generalplot.R")

plot_coverage <- function(filename, outputname){
    # Plot coverage histogram and coverage itself
    cov <- read_tsv(filename,col_names=FALSE, col_types = 'cii')
    names(cov) <- c("name","index","cov")
    g <- cov %>% ggplot(aes(x = cov)) + geom_histogram(bins = 100) +
        labs(x = "coverage", y = "count", title = filename)
    generalplot(g,paste0(outputname,"histogram"))    
    count <- cov %>% count() %>% pull(n)
    bins <- 1000
    binwidth <- ceiling(count / bins) 
    data <-  cov %>% mutate(bin_id = rep(1:bins,each=binwidth)[1:count]) %>% nest(-bin_id) %>%
        mutate(data = map(data,function(x) x %>% summarize(mean = mean(cov),name = name[1]))) %>%
        unnest()
    g <- data %>% mutate(position = bin_id * binwidth) %>% 
        ggplot(aes(x = position, y = mean)) + geom_line() 
    generalplot(g,paste0(outputname,"positionwize"))
}


sequel <- "/grid/ban-m/arabidopsis_thaliana/sequel/coverage/"
nanopore <- "/grid/ban-m/arabidopsis_thaliana/nanopore/coverage/"

coverages <- list(list("coverage_mito_MAPQ0.wig","coverage_mito_mapq0"),
                  list("coverage_mito_MAPQ10.wig","coverage_mito_mapq10"),
                  list("coverage_1_MAPQ0.wig","coverage_1_mapq0"),
                  list("coverage_1_MAPQ10.wig","coverage_1_mapq10"))

coverages %>% sapply(FUN=function(ls){
    plot_coverage(paste0(sequel,"sequel_minialign_",ls[[1]]),
                  paste0("sequel_minialign_",ls[[2]]))})

coverages %>% sapply(FUN=function(ls){
    plot_coverage(paste0(nanopore,"nanopore_minialign_",ls[[1]]),
                  paste0("nanopore_minialign_",ls[[2]]))})



## ------------------Pull the extremely low coverage --------------------------------

sequeldata <- read_tsv("/grid/ban-m/arabidopsis_thaliana/sequel/coverage/sequel_minialign_coverage_mito.wig",
                       col_names = FALSE,
                       col_types = 'cii') %>% arrange(X2)

threshold <- sequeldata %>% pull(X3) %>% quantile(.05)

low_coverage <- sequeldata %>% filter(X3 < threshold)
index <- low_coverage %>% pull(X2)
index <- index + 1
index <- c(0,index[1:dim(low_coverage)[1]-1])
group <- cumsum(index)

low_coverage_region <- low_coverage %>%
    mutate(class = X2 - index) %>%
    mutate(group = cumsum(class)) %>%
    nest(-group) %>%
    mutate(size = map(data,function(x) x %>% count() %>% pull(n))) %>%
    unnest(size) %>% filter(1000 < size)

result <- low_coverage_region %>%
    mutate(data = map(data,function(x) x  %>% summarize(start = min(X2),
                                                        end = max(X2)))) %>%
    unnest()
write_tsv(result,"./data/sequel_minialign_low_coverage_region.tsv")
