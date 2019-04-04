library("tidyverse")
library("stringi")
args <- c("./result/sequel_positionwise_coverage_short.tsv",
          "./result/sequel_minimap2_assignment.tsv")

positionwise_coverage <- readLines(args[1])
assignment <- read_tsv(args[2],col_name = FALSE, col_types = "cc") %>%
    rename(id = X1, type = X2)

assign_table <- assignment %>% mutate(chrtype = ifelse(stri_detect_charclass(type,"[0-9]"),"genome",type)) %>%
    pull(chrtype)
names(assign_table) <- assignment$id

to_vectors <- function(input){
    strsplit(x = input, split = ',') %>% lapply(function(ls){
        name <- ls[1]
        value <- as.numeric(ls[-1])
        list(name = name,
             value = value[!is.na(value)],
             type = assign_table[name]
             )
    })
}

positionwise_parsed <- to_vectors(positionwise_coverage)

cal_mean <- function(ls){
    list(name = ls$name,
         mean = mean(ls$value))
}

to_tibble <- function(ls){
    tibble(coverage = ls$value,
           index = 1:length(ls$value))
}

test <-  positionwise_parsed[[1]] %>% to_tibble()

thr <- (test %>% pull(coverage) %>% mean()) * 0.03
mean <-  test %>% pull(coverage) %>% mean()
sd <-  test %>% pull(coverage) %>% sd()
thr <- mean - 3 * sd 

len <- test %>% pull(coverage) %>% length() 
start <- floor(len*0.03)
end <- floor(len*0.97)
g1 <- qplot(x = index, y = coverage, data = test)
g2 <- qplot(x = index, y = coverage, data = test[start:end,])
