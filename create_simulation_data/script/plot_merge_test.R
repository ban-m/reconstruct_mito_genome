library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

dataset <- read_tsv("./result/merge_test_diff.txt",col_names=FALSE) %>%
    rename(seed = X1, lkgain = X2, dist = X3)
g <- dataset %>% filter(dist < 5) %>% 
    ggplot() + geom_point(mapping = aes(x = dist, y = lkgain))


generalplot(g,"merge_criteria")
