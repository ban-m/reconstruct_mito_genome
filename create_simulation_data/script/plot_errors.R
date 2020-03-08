library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
## generalplot <- function(g,name){
##     cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
##                     plot = g + cowplot::theme_cowplot(font_size=12),
##                     dpi=350,width = 178,height = 86,units="mm")
##     cowplot::ggsave(filename = paste0("./png/",name,".png"),
##                     plot = g + cowplot::theme_cowplot(font_size=12),
##                     dpi = 350,width = 178,height = 86,units="mm")
## }

generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}


args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
## dataset <- read_tsv("./result/last_decompose_gibbs_test.txt", col_names=FALSE)
## dataset <- read_tsv("./result/last_decompose_gibbs_150.tsv", col_names=FALSE)
## dataset <- read_tsv("./result/last_decompose_poa.tsv", col_names=FALSE)
## dataset <- read_tsv("./result/last_decompose_dbg.tsv", col_names=FALSE)
dataset <- read_tsv(filename, col_names=FALSE)

accs <- dataset %>% mutate(acc1 = X4 + X7, acc2 = X5 + X6, acc = pmax(acc1, acc2)) %>%
    mutate(acc = acc / X9) %>% select(X2,X3,acc1,acc2, acc, X8,X9,X10)


mean_accs <-  accs %>% mutate(Coverage = X2 + X3) %>%
    select(Coverage, X10, acc) %>% nest(-X10, -Coverage) %>%
    mutate(data = map(data, ~summarize(.,Accuracy = median(acc)))) %>%
    unnest() 

g <- mean_accs %>%
    rename(DivergenceRate = X10) %>%
    mutate(DivergenceRate = 100 * DivergenceRate) %>% 
    ggplot() +
    geom_raster(mapping = aes(x = Coverage, y = DivergenceRate, fill = Accuracy))  +
    labs(y = "Divergence Rate(%)") +
    scale_fill_continuous(limits =c(0.5,1))

generalplot(g,paste0(outname))

