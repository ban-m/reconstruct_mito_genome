library(tidyverse)
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi=350,width = 178,height = 86,units="mm")
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi = 350,width = 178,height = 86,units="mm")
}

dataset <- read_tsv("./result/entropies.tsv",col_names=FALSE)

g <- dataset %>% nest(-X1) %>% mutate(Replicate = row_number()) %>%
    mutate(data = map(data, ~ mutate(., Iteration = row_number()))) %>% unnest() %>%
    ggplot() + geom_line(mapping = aes(x = Iteration, y = X2, color = factor(Replicate))) +
    labs(y = "Total Entropy")
generalplot(g,"entropy_trajectry")


g <- dataset %>% nest(-X1) %>% mutate(Replicate = row_number()) %>%
    mutate(data = map(data, ~ mutate(., Iteration = row_number()))) %>% unnest() %>%
    ggplot() +
    geom_line(mapping = aes(x = Iteration, y = pmax(1-X9,X9), color = factor(Replicate))) +
    labs(y = "Accuracy")
generalplot(g,"accuracy_trajectry")
