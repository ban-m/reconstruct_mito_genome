library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

filename <- "./result/variable_chain_length.tsv"


dataset <- read_tsv(filename)

g <- dataset %>% ggplot() + geom_point(aes(x = Dist/Length, y = HMM, color = Length))
generalplot(g,"vlc_point")

dataset %>% select(-Dist, Coverage) %>% nest(-Length)  %>%
    mutate(data= map(data,function(x) x %>% summarize(mean = mean(HMM))))


g <- dataset %>% ggplot(mapping = aes(x = Length, y = HMM))  + geom_point()
generalplot(g,"vlc_point_2")
g <- dataset %>% ggplot(mapping = aes(x = Length, y = HMM))  + geom_smooth()
generalplot(g,"vlc_smooth_2")
