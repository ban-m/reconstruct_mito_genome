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
    mutate(Replicate = factor(row_number())) %>%
    mutate(data = map(data, function(x) x %>% mutate(Iteration = 1:dim(x)[1]))) %>% 
    unnest() %>%
    rename(LogLikelihood = X2) 

g <- dataset %>% ggplot() + geom_line(aes(x = Iteration, y = LogLikelihood, color = Replicate))
generalplot(g, outputname)
