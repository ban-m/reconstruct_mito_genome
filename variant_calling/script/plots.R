library("tidyverse")
args <- commandArgs(trailingOnly=TRUE)
loadNamespace("cowplot")
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}
filename <- args[1]
outputname <- args[2]
data <- read_tsv(filename)

g <-data %>%  filter(MAF < 0.2) %>%
    ggplot() +
    geom_histogram(mapping = aes(x = MAF, y = ..density..), bins=100) + facet_grid(Type ~ .)
    stat_function(fun = function(x){dnorm(x,mean=mean(data$MAF),sd=sd(data$MAF))})
ggsave(filename = str_c("./pics/", outputname, "_histogram.png"), g)

g <-data %>%  filter(MAF<0.5) %>% ggplot() +
    geom_point(mapping = aes(x = Pos, y = MAF, color = Type), size = 0.2) + facet_grid(Type ~ .) + ylim(c(0,0.5))
generalplot(g,str_c(outputname,"_point"))

g <-data %>%  filter(MAF<0.5) %>% ggplot() +
    geom_point(mapping = aes(x = Pos, y = MAF, color = Type), size = 0.2) + ylim(c(0,0.5))
generalplot(g, str_c(outputname,"_point_merged"))


data <- data %>% filter(Type == "Sub")
filtered <- data %>%
    filter((4956 < Pos & Pos < 37653) |
           (109386 < Pos & Pos < 144786) |
           (341256 < Pos & Pos < 367807) )
(37653 - 4956 + 144786 - 109386 + 367807 - 341256)/367807
print("Threshold")
print(mean(data$MAF) + 2*sd(data$MAF))
big_maf <- data %>% filter(MAF > mean(data$MAF) + 2*sd(data$MAF))
big_maf_count <- big_maf %>% count() %>% pull(n)
total <- data %>% count() %>% pull(n)
print("Big MAF Pisition(Exclude putative nucleic region.)")
print(big_maf_count)
print("Total number of position investigated.")
print(total)
print("Fraction of BigMaf")
print(big_maf_count/total * 100)
