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

dataset <- read_tsv("./result/maximum_rank_min.tsv")
mindata <- dataset %>% nest(-Coverage, -Length)  %>% mutate(data = map(data,~slice(., which.min(Min)))) %>% unnest()
g <- dataset %>% filter(Length == 150) %>% ggplot() + geom_point(aes(x = Coverage, y = Min))
g2 <- mindata %>% filter(Length == 150) %>% ggplot() + geom_point(aes(x = Coverage, y = Min))

result <- lm(Min~Coverage,mindata %>% filter(Length == 150 & Coverage < 50) )

dataset_rowcov <- read_tsv("./result/maximum_rank_min_lowcov.tsv")
mindata_rowcov <- dataset_rowcov %>%
    nest(-Coverage)  %>%
    mutate(data = map(data,~slice(., which.min(Min)))) %>%
    unnest()
result_rowcov <- lm(Min~Coverage, data= mindata_rowcov %>% filter(Coverage < 21))
