library(tidyverse)
data <- read_tsv("./result/real_dataset.tsv", col_names=FALSE)
data <- data %>% select(-X4) %>% rename(Match=X1, Ins=X2,  Del=X3)
data %>% summarize(Match_m= mean(Match), Match_sd = sd(Match),
                   Del_m = mean(Del), Del_sd = sd(Del),
                   Ins_m = mean(Ins), Ins_sd = sd(Ins),
                   )
