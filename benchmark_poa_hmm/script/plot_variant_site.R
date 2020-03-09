library("tidyverse")
dataset <- read_tsv("./result/variant_calling.tsv", col_names=FALSE)
g <- dataset %>% mutate(X3 = factor(X3)) %>%
    rename(NumberOfVatiants= X3, Position = X2, Weight = X4) %>% 
    ggplot() +
    geom_point(mapping = aes(x = Position, y =Weight, color = NumberOfVatiants), size =3) +
    facet_grid(X5 ~ X1)
ggsave("./png/variant_call.png", g)
ggsave("./pdf/variant_call.png", g)

