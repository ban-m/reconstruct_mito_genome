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
dataset <- read_tsv(filename)

ds <- dataset %>% gather(key = Type, value = Accuracy, -Test, -Dist, -Training)

g <- ds %>%
    ggplot() +
    geom_point(mapping = aes(x = Test, y = Accuracy, color = Type), alpha = 0.4) + facet_wrap(Dist~.)
generalplot(g, paste0(outputname,"_point"))

g <- ds %>%
    ggplot() +
    geom_histogram(mapping = aes(x =Accuracy, y = ..density.., fill = Type), position = "identity",alpha = 0.4) 
generalplot(g, paste0(outputname,"_hist_no_dist"))

g <- ds %>%
    ggplot() +
    geom_histogram(mapping = aes(x =Accuracy, y = ..density.., fill = Type), position = "identity",alpha = 0.4)  +
    facet_wrap(Dist~.)
generalplot(g, paste0(outputname,"_hist"))


    

g <- ds %>%
    ggplot() +
    geom_smooth(mapping = aes(x = Test, y = Accuracy, color = Type)) + facet_wrap(Dist~.)
generalplot(g, paste0(outputname,"_smooth"))

g <- dataset %>% mutate(Improvement = EM- Naive) %>%
    ggplot() +
    geom_violin( mapping = aes (x = factor(Test), y = Improvement))
generalplot(g, paste0(outputname,"_violin"))

g <- dataset %>% mutate(Improvement = EM- Naive) %>%
    filter(Test > 50) %>% 
    ggplot() +
    geom_violin(mapping = aes (x = factor(Training), y = Improvement))
generalplot(g, paste0(outputname,"_violin_Training"))

matrix <- dataset %>%
    mutate(Improvement = EM- Naive) %>%
    select(Improvement, Training, Test) %>%
    nest(Improvement) %>%
    mutate(data = map(data, function(x) x %>% summarize(mean = mean(Improvement)))) %>%
    unnest()
g <- matrix %>%
    ggplot() +
    geom_raster(mapping = aes(x = Training, y = Test, fill = mean))
generalplot(g, paste0(outputname,"_raster"))

