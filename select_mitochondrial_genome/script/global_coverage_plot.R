library("tidyverse")
library("stringi")
loadNamespace("cowplot")
args <- commandArgs(trailingOnly = TRUE)

### ------ Load data -----

args <- c("./result/sequel_average_coverage_trimed3.tsv",
          "./result/sequel_minimap2_assignment.tsv")

args <- c("./result/sequel_average_coverage.tsv",
          "./result/sequel_minimap2_assignment.tsv")

data <- read_tsv(args[1],col_names = FALSE) %>%
    rename(id = X1, coverage = X2, variance = X3, length = X4)
annot <- read_tsv(args[2],col_names = FALSE,col_types = c("cc")) %>%
    rename(id = X1, type = X2)
data <- full_join(data,annot) %>% filter(!is.na(coverage))  %>%
    mutate(chrtype = ifelse(stri_detect_charclass(type,"[0-9]"),"genome",type)) 


g <- data %>% filter(coverage < 75) %>%
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(mapping=aes(x = coverage,fill = chrtype),
                   position = "identity",
                   alpha = 0.6,
                   bins = 60) +
    cowplot::theme_cowplot() 
cowplot::ggsave(filename = "./pdf/sequel_coverage_UB_75.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_UB_75.png",plot=g)

g <- data %>% filter(coverage < 150) %>%
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(mapping=aes(x = coverage,fill = chrtype),
                   position = "identity",
                   alpha = 0.6,
                   bins = 60) +
    cowplot::theme_cowplot() 
cowplot::ggsave(filename = "./pdf/sequel_coverage_UB_150.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_UB_150.png",plot=g)


g <- data %>% 
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(mapping=aes(x = coverage,fill = chrtype),
                   position = "identity",
                   alpha = 0.6,
                   bins = 60) +
    cowplot::theme_cowplot() 
cowplot::ggsave(filename = "./pdf/sequel_coverage.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage.png",plot=g)


g <- data %>% filter(coverage > 50) %>% 
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(mapping=aes(x = coverage,fill = chrtype),
                   position = "identity",
                   alpha = 0.6,
                   bins = 60) +
    cowplot::theme_cowplot() 
cowplot::ggsave(filename = "./pdf/sequel_coverage_LB_50.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_LB_50.png",plot=g)

g <- data %>% filter(coverage > 50 & coverage < 1000) %>% 
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(mapping=aes(x = coverage,fill = chrtype),
                   position = "identity",
                   alpha = 0.6,
                   bins = 60) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_LB_100_UP_1000.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_LB_100_UP_1000.png",plot=g)


g <- data %>% filter(coverage > 100 & coverage < 1000) %>% 
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(bins = 60) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_LB_100_UP_1000_wc.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_LB_100_UP_1000_wc.png",plot=g)

g <- data %>% 
    ggplot(mapping = aes(x = coverage)) + 
    geom_histogram(bins = 60) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_wc.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_wc.png",plot=g)



g <- data %>%
    ggplot(mapping = aes(x = coverage, y = variance, color = chrtype)) +
    geom_point(alpha = 0.4) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_scatter.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_scatter.png",plot=g)

g <- data %>% filter(50 < coverage & coverage < 1000 & variance < 500) %>% 
    ggplot(mapping = aes(x = coverage, y = variance, color = chrtype)) +
    geom_point(alpha = 0.2) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_scatter_zoom1.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_scatter_zoom1.png",plot=g)



g <- data %>% filter(250 < coverage & coverage < 750 & variance > 20 & variance < 200) %>%
    ggplot(mapping = aes(x = coverage, y = variance, color = chrtype)) +
    geom_point(alpha = 0.2) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_scatter_zoom2.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_scatter_zoom2.png",plot=g)

g <- data %>%
    filter(chrtype == "mitochondria") %>% 
    filter(250 < coverage & coverage < 750 & variance > 20 & variance < 200) %>%
    ggplot(mapping = aes(x = coverage, y = variance, color = chrtype)) +
    geom_point(alpha = 0.2) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/sequel_coverage_scatter_zoom3.pdf",plot=g)
cowplot::ggsave(filename = "./png/sequel_coverage_scatter_zoom3.png",plot=g)

