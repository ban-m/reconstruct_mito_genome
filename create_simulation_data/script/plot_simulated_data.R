library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outputname <- args[2]

dataset <- read_tsv(filename) %>% 
    gather(key = Type, value = Accuracy, -Dist, -Coverage, -Length) %>%
    mutate(ErrorRate = Dist/Length * 100)

g <- dataset %>%
    ggplot() +
    geom_point(mapping = aes(y = Accuracy, x = ErrorRate, color = Type), alpha=0.4) +
    facet_wrap(.~Coverage)
generalplot(g, paste0(outputname,"_point"))


g <- dataset %>%
    ggplot() +
    geom_smooth(mapping = aes(y = Accuracy, x = ErrorRate, color = Type), alpha=0.4) +
    facet_wrap(.~Coverage)
generalplot(g, paste0(outputname,"_smooth"))

temp <- dataset %>% select(-ErrorRate, -Length)

temp %>%
    nest(Accuracy) %>%
    mutate(data = map(data, function(x) x %>% summarize(mean = mean(Accuracy)))) %>%
    unnest() %>%
    nest(mean,Coverage, Dist) %>%
    mutate(data = map(data, function(x) x %>% summarize(area = sum(mean)))) %>%
    unnest() %>%
    spread(key = Type, value = area)
