library("tidyverse")
args<- commandArgs(trailingOnly = TRUE)
                                        # args<- c("./result/tiling_subunit_count.tsv")

countdata <- read_tsv(args[1])

g <- countdata %>% mutate(type = as.factor(paste0(ctg,":",unit))) %>% 
    ggplot(mapping= aes(x = subunit,y=count)) + geom_line() + facet_wrap(type ~ . )
ggsave(filename = "./png/unit_count.png",plot=g)
