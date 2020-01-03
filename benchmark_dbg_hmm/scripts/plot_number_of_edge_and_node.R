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
dataset <- read_tsv("./result/number_of_edge_and_node.tsv")
g <- dataset %>% filter(Coverage < 200) %>% 
    gather(key = Type, value = Count, -Coverage, -Weight, -Method) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Count, color = Type)) +
    facet_wrap(.~Method)
generalplot(g,"number_of_edge_and_node")

g <- dataset %>% filter(Coverage < 30) %>% 
    gather(key = Type, value = Count, -Coverage, -Weight, -Method) %>%
    ggplot() + geom_point(aes(x = Coverage, y = Count, color = Type)) +
    facet_wrap(.~Method)
generalplot(g,"number_of_edge_and_node_30")
