library("tidyverse")
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi=350,width = 178,height = 86,units="mm")
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi = 350,width = 178,height = 86,units="mm")
}

## generalplot <- function(g,name){
##     cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
##                     plot = g + cowplot::theme_cowplot())
##     cowplot::ggsave(filename = paste0("./png/",name,".png"),
##                     plot = g + cowplot::theme_cowplot())
## }

dataset <- read_tsv("./result/neutral_prediction.tsv")
g <- dataset %>% ggplot() + geom_histogram(mapping = aes(x = Entropy),bins = 100) + facet_grid( Skew ~ Type)
x <- dataset %>% filter(Type == "HMM" & Skew == "Skew") %>% pull(Entropy) %>% summary()
print(x)
count <- dataset %>% filter(Type == "HMM" & Skew == "Neutral") %>% filter(Entropy > x["3rd Qu."]) %>% count()
print("total")
print(dataset %>% filter(Type == "HMM" & Skew == "Neutral") %>% count())
print("Below, more than 3rd quantile")
print(count)

generalplot(g, "entropies")
