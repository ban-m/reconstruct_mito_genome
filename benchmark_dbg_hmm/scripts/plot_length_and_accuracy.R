library("tidyverse")
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
dataset <- read_tsv(args[1])
dataset <- dataset %>% gather(key = type, value = accuracy, -Length, -Dist)


g <- dataset %>% ggplot() + geom_smooth(mapping = aes( x = Length, y = accuracy, color = type), se = TRUE) + facet_wrap( . ~ Dist)
generalplot(g, paste0("gam_",outname))


g <- dataset %>%
    ggplot() +
    geom_smooth(mapping = aes( x = Length, y = accuracy, linetype = type), se = TRUE, color = "black") + facet_wrap( . ~ Dist)
generalplot(g, paste0("gam_",outname,"_bw"))


g <- dataset %>% ggplot() + geom_point(mapping = aes( x = Length, y = accuracy, color = type), alpha = 0.09) + facet_wrap( . ~ Dist)
generalplot(g, paste0("point_",outname))

summary <- dataset %>% nest(accuracy) %>%
    mutate(data = map(data,
                      function(x) x  %>% summarize(mean = mean(accuracy),
                                                   sd = sd(accuracy)))) %>%
    unnest()
write_tsv(summary, path = paste0("./result/",outname,"_summary.tsv"))
