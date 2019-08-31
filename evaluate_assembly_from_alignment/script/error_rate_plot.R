library("tidyverse")
error_rate <- read_tsv("./result/mitochondria/error_profile.dat",col_names=FALSE) %>%
    rename(index = X1,
           subst = X2,
           ins = X3,
           del = X4,
           depth = X5)
len <- length(error_rate$index)
formatted <- tibble(index = rep(error_rate$index,4),
                    values = c(error_rate$subst,error_rate$ins,error_rate$del,error_rate$depth),
                    type = c(rep("Subst",len),rep("Ins",len),rep("Del",len),rep("Depth",len)))
g <-formatted %>%  ggplot(mapping = aes(x = index, y = values, color = type)) + geom_line()
