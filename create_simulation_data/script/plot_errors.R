library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi=600,width = 178,height = 86,units="mm")
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi = 600,width = 178,height = 86,units="mm")
}

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outname <- args[2]
# dataset <- read_tsv("./result/last_decompose_num_poa.tsv", col_names=FALSE)
dataset <- read_tsv(filename, col_names=FALSE)

accs <- dataset %>% mutate(acc1 = X4 + X7, acc2 = X5 + X6, acc = pmax(acc1, acc2)) %>%
    mutate(acc = acc / X9) %>% select(X2,X3,acc1,acc2, acc, X8,X9,X10)


mean_accs <-  accs %>% mutate(Coverage = (X2 + X3)/2) %>%
    select(Coverage, X10, acc) %>% nest(-X10, -Coverage) %>%
    mutate(data = map(data, ~summarize(.,Accuracy = mean(acc)))) %>%
    unnest() 

g <- mean_accs %>%
    rename(VariantSite = X10) %>%
    ggplot() +
    geom_raster(mapping = aes(x = Coverage, y = VariantSite, fill = Accuracy))  +
    labs(y = "Number of Variant Sites") +
    scale_fill_continuous(limits =c(0.5,1))
generalplot(g,paste0(outname))

g <- accs %>%
    filter(X10 <= 4) %>%
    ggplot() + geom_smooth(mapping = aes(x = X9, y = acc, color = factor(X10))) + 
    labs(x = "Total Coverage", y = "Accuracy")
generalplot(g,paste0(outname,"_smooth"))

g <- accs %>% 
    nest(-X8) %>%
    mutate(data = map(data,~summarize(.,Accuracy = mean(acc)))) %>%
    unnest() %>% 
    ggplot() + geom_point(mapping = aes(x = X8, y = Accuracy))
generalplot(g,paste0(outname,"_point"))

g <- accs %>%
    filter(X10 <= 4) %>%
    ggplot() + geom_boxplot(mapping = aes(x = factor(X9), y = acc, color = factor(X10))) +
    labs(x = "Total Coverage", y = "Accuracy", color = "Number of Variants\nout of 9K bp")

generalplot(g,paste0(outname,"_boxplot"))

g <- accs %>% filter(X9 %in% c(80, 110, 140)) %>%
    ggplot() + geom_boxplot(mapping = aes(x = factor(X9), y = acc, color = factor(X10))) +
    labs(x = "Total Coverage", y = "Accuracy", color = "Number of Variants\nout of 9K bp") +
    coord_flip()
generalplot(g,paste0(outname,"_boxplot_flip"))
