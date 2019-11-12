library("tidyverse")
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

dataset <- read_tsv("./result/neutral_prediction.tsv")
g <- dataset %>% ggplot() + geom_histogram(mapping = aes(x = Entropy),bins = 100) + facet_grid( Skew ~ Type)
generalplot(g, "entropies")
