library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

dataset <- read_tsv("./result/simulated_data.tsv") %>%
    gather(key = Type, value = Accuracy, -Dist, -Coverage, -Length) %>%
    mutate(ErrorRate = Dist/Length * 100)

g <- dataset %>%
    ggplot() +
    geom_point(mapping = aes(y = Accuracy, x = ErrorRate, color = Type), alpha=0.4) +
    facet_wrap(.~Coverage)
generalplot(g, "simulated_data_point")


g <- dataset %>%
    ggplot() +
    geom_smooth(mapping = aes(y = Accuracy, x = ErrorRate, color = Type), alpha=0.4) +
    facet_wrap(.~Coverage)
generalplot(g, "simulated_data_smooth")
