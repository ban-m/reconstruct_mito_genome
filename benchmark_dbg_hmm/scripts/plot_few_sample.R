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
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outputname <- args[2]
dataset <- read_tsv(filename) %>%
    select(-HMM) %>% 
    gather(key = Type, value = Accuracy, -Dist, -Coverage) 

upperbound <- dataset %>% select(Dist, Coverage) %>%
    mutate(UpperBound = pmax(1-0.14**Dist,0.5))


summary <- dataset %>%
    filter(Type == "WHMM") %>%
    nest(-Dist, -Coverage) %>%
    mutate(data = map(data,function(x)summarize(x, mean = mean(Accuracy)))) %>%
    unnest() %>% arrange(Coverage, Dist) 

print(summary, n = Inf)


g <- dataset %>%
    ggplot() +
    geom_point(mapping = aes(x = Coverage, y = Accuracy), alpha = 0.5, size = 1, stroke = 0) + 
    geom_line(mapping = aes(x = Coverage, y = UpperBound), data = upperbound)  + 
    facet_grid(Type~Dist)

generalplot(g,paste0(outputname,"_point"))

g <- dataset %>%
    ggplot() +
    geom_smooth(mapping = aes(x = Coverage, y = Accuracy, color = Type)) +
    geom_line(mapping = aes(x = Coverage, y = UpperBound), data = upperbound)  + 
    facet_wrap(.~Dist)
generalplot(g,paste0(outputname,"_smooth"))

g <- dataset %>%
    filter(Dist > 0 & Dist < 4) %>% 
    ggplot() +
    geom_bin2d(mapping = aes(x = Coverage, y = Accuracy)) +
    scale_fill_gradient(low = "white", high = "black") + 
    facet_grid(Type~Dist)


dataset <- read_tsv(filename)
g <-  dataset %>%
    mutate(ImpHMM = WHMM-HMM, ImpAln = WHMM-Aln) %>%
    filter(Dist > 0) %>%
    select(-WHMM, -HMM, -Aln) %>% 
    gather(key = Type, value = Improvement,-Dist,-Coverage ) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Improvement, color = Type))  + facet_wrap(.~Dist)

g <-  dataset %>%
    mutate(Improvement = WHMM-Aln) %>%
    filter(Dist > 0) %>%
    select(-WHMM, -HMM, -Aln) %>% 
    gather(key = Type, value = Improvement,-Dist,-Coverage ) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Improvement))  + facet_wrap(.~Dist)
generalplot(g,paste0(outputname,"_improvement_point"))

temp <- dataset %>% mutate(Improvement = WHMM - Aln) %>%
    select(Improvement, Coverage) %>% 
    nest(-Coverage) %>%
    mutate(data = map(data,function(x) x %>% summarize(mean = mean(Improvement),
                                                       sd = sd(Improvement),
                                                       max = max(Improvement)))) %>% 
    unnest()
