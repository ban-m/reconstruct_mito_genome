library("tidyverse")
dataset <- read_tsv("./lks", col_names = FALSE) %>%
    rename(Read=X1,Type = X2, Index=X3, LK = X4)



test_pair <- function(i,j) {
    temp <- dataset %>% filter(Index %in% c(i,j) & Type == 0) %>%
        mutate(Index = ifelse(Index == i, "axis_i", "axis_j")) %>%
        spread(key = Index,value = LK )
    summary(lm(axis_i ~ axis_j, data =temp))$coefficients[2,4]
}


lk_test <- function(df){
    wilcox.test(formula = LK ~ Type, data = df)
}

result <- dataset %>% nest(-Index) %>%
    mutate(data = map(data, lk_test)) %>%
    mutate(data = map(data, function(c) c$p.value)) %>% unnest()
    
result <- tibble::new_tibble(data.frame(t(combn(x = 0:39, m=2)))) %>%
    mutate(pvalue = mapply(test_pair, X1,X2))
