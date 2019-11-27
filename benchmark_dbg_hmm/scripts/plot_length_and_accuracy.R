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
dataset <- read_tsv(filename) %>%
    gather(key = type, value = accuracy, -Length, -Dist)

upperbound <- dataset %>%
    select(Dist, Length) %>% mutate(UpperBound = 1-0.14**Dist) 

g <- dataset %>% ggplot() +
    geom_smooth(mapping = aes( x = Length, y = accuracy, color = type), se = TRUE) +
    geom_line(mapping = aes(x = Length, y = UpperBound), data = upperbound) +
    facet_wrap( . ~ Dist) 
generalplot(g, paste0("gam_",outname))


g <- dataset %>%
    ggplot() +
    geom_smooth(mapping = aes( x = Length, y = accuracy, linetype = type), se = TRUE, color = "black") +
    geom_line(mapping = aes(x = Length, y = UpperBound), data = upperbound, linetype = "dotted") +
    facet_wrap( . ~ Dist)
generalplot(g, paste0("gam_",outname,"_bw"))


g <- dataset %>% filter(Type == "Proposed") %>%
    mutate(ErrorRate = Dist/Length) %>% 
    ggplot() +
    geom_line(mapping = aes(x = Length, y = UpperBound), data = upperbound) +
    geom_smooth(mapping = aes( x = ErrorRate, y = accuracy, color = factor(Length)), alpha = 0.09)
generalplot(g, paste0("gam_errorrate_",outname))

g <- dataset %>% filter(Type == "Proposed") %>%
    mutate(ErrorRate = Dist/Length) %>% 
    ggplot() + geom_point(mapping = aes( x = ErrorRate, y = accuracy, color = Length), alpha = 0.09)
generalplot(g, paste0("gam_errorrate_",outname))


g <- dataset %>% ggplot() + geom_point(mapping = aes( x = Length, y = accuracy, color = type), alpha = 0.09) + facet_wrap( . ~ Dist)
generalplot(g, paste0("point_",outname))

summary <- dataset %>% nest(accuracy) %>%
    mutate(data = map(data,
                      function(x) x  %>% summarize(mean = mean(accuracy),
                                                   sd = sd(accuracy)))) %>%
    unnest()
write_tsv(summary, path = paste0("./result/",outname,"_summary.tsv"))
