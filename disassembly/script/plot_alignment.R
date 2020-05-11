loadNamespace("cowplot")
library(tidyverse)
## ==== Setting for publishing figure ====== 
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pics/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi=350,width = 178,height = 86,units="mm")
    cowplot::ggsave(filename = paste0("./pics/",name,".png"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi = 350,width = 178,height = 86,units="mm")
}


args <- commandArgs(trailingOnly = TRUE)
dataset <- read_tsv(args[1])
g <- dataset %>% ggplot() +
    geom_tile(mapping = aes(x = Query, y = Target, fill = CoverRate)) +
    labs(x = "Query Accession", y = "Target Accession", fill = "Cover Rate")
generalplot(g, args[2])
