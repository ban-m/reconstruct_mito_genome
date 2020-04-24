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
dataset <- read_tsv(args[1], col_names=FALSE)
g <- dataset %>% ggplot() + geom_tile(mapping = aes(x = X1, y = X3, fill = X2)) +
    labs(x = "names of PacBio contig", y = "Accession ID")
generalplot(g, args[2])
