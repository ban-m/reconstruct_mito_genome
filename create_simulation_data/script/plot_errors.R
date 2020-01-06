library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
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



dataset <- read_tsv("./result/last_decompose_gibbs_ve.tsv", col_names=FALSE)

accs <- dataset %>% mutate(acc1 = X4 + X7, acc2 = X5 + X6, acc = pmax(acc1, acc2)) %>%
    mutate(acc = acc / X9) %>% select(X2,X3,acc1,acc2, acc, X8,X9)


mean_accs <-  accs %>% nest(-X9,-X8) %>%
    mutate(data = map(data, function(x)summarize(x,mean = mean(acc)))) %>%
    unnest() 

tiled_accs <- accs %>%
    mutate(Coverage = X9 %/% 10 * 10, NumOfError = (X8 %/% 5 * 5)) %>%
    nest(-NumOfError, -Coverage) %>%
    mutate(data = map(data, function(x)summarize(x,mean = mean(acc)))) %>%
    unnest()

g <- tiled_accs %>% ggplot() + geom_raster(mapping = aes(x = Coverage, y = NumOfError, fill = mean))
generalplot(g,"150_40_varying_error_rate")

g <- accs %>% nest(-X9,-X8) %>%
    mutate(data = map(data, function(x)summarize(x,mean = mean(acc)))) %>%
    unnest() %>%
    ggplot() + geom_raster(mapping = aes(x = X8, y = X9, fill = mean), interpolate = TRUE)

g <- accs %>%
    mutate(ErrorRate = X8 / 150 / 40) %>% filter(X9 > 100) %>% 
    ggplot() + geom_hex(mapping = aes(x = ErrorRate, y = acc))
