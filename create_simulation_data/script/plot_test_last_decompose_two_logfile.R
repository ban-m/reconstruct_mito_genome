library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}


args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
outputname <- args[2]

dataset <- read_tsv(filename,col_names = FALSE) %>% 
    nest(-X1) %>% select(-c(1)) %>% 
    mutate(class = factor(row_number())) %>%
    mutate(data = map(data, function(x) x %>% mutate(iter = 1:dim(x)[1]))) %>% 
    unnest()

g <- dataset %>% ggplot() + geom_line(aes(x = iter, y = X2, color = class))
generalplot(g, outputname)
