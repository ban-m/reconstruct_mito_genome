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
outputname <- args[2]
dataset <- read_tsv(filename) %>% gather(key = Type, value = Accuracy, -Dist, -Coverage)



upperbound <- dataset %>% select(Dist, Coverage) %>%
    mutate(UpperBound = pmax(1-0.14**Dist,0.5))

g <- dataset %>%
    ggplot() +
    geom_point(mapping = aes(x = Coverage, y = Accuracy, color = Type), alpha = 0.5) +
    geom_line(mapping = aes(x = Coverage, y = UpperBound), data = upperbound)  + 
    facet_wrap(.~Dist)
generalplot(g,paste0(outputname,"_point"))

g <- dataset %>%
    ggplot() +
    geom_smooth(mapping = aes(x = Coverage, y = Accuracy, color = Type)) +
    geom_line(mapping = aes(x = Coverage, y = UpperBound), data = upperbound)  + 
    facet_wrap(.~Dist)
generalplot(g,paste0(outputname,"_smooth"))
