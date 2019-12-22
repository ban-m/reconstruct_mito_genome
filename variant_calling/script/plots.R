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
data <- read_tsv(filename, col_names =FALSE) %>% rename(Position=X1, Freq=X2)

g <-data %>%  ggplot() +
    geom_histogram(mapping = aes(x = Freq, y = ..density..), bins=100) +
    stat_function(fun = function(x){dnorm(x,mean=mean(data$Freq),sd=sd(data$Freq))})
ggsave(filename = str_c("./pics/", outputname, "_histogram.png"), g)

g <-data %>%  ggplot() +
    geom_point(mapping = aes(x = Position, y = Freq))
ggsave(filename = str_c("./pics/",outputname,"_point.png"), g)

filtered <- data %>%
    filter((4956 < Position & Position < 37653) |
           (109386 < Position & Position < 144786) |
           (341256 < Position & Position < 367807) )
(37653 - 4956 + 144786 - 109386 + 367807 - 341256)/367807
print("Threshold")
print(mean(data$Freq) + 3*sd(data$Freq))
big_maf <- data %>% filter(Freq > mean(data$Freq) + 3*sd(data$Freq))
big_maf_count <- big_maf %>% count() %>% pull(n)
total <- data %>% count() %>% pull(n)
print("Big MAF Pisition(Exclude putative nucleic region.)")
print(big_maf_count)
print("Total number of position investigated.")
print(total)
print("Fraction of BigMaf")
print(big_maf_count/total * 100)
