library("tidyverse")
library("stringi")
loadNamespace("cowplot")
args <- commandArgs(trailingOnly = TRUE)

### ------ Load data -----
## args <- c("./result/sequel_positionwise_coverage.tsv",
##           "./result/sequel_minimap2_assignment.tsv")

positionwise_coverage <- readLines(args[1])
assignment <- read_tsv(args[2],col_name = FALSE, col_types = "cc") %>% rename(id = X1, type = X2)

assign_table <- assignment %>% mutate(chrtype = ifelse(stri_detect_charclass(type,"[0-9]"),"genome",type)) %>%
    pull(chrtype)
names(assign_table) <- assignment$id

to_vectors <- function(input){
    strsplit(x = input, split = ',') %>% lapply(function(ls){
        name <- ls[1]
        value <- as.numeric(ls[-1])
        list(name = name,
             value = value[!is.na(value)],
             type = assign_table[name]
             )
    })
}

positionwise_parsed <- to_vectors(positionwise_coverage)

cal_mean <- function(ls){
    list(name = ls$name,
         mean = mean(ls$value))
}


chrname_to_color <- function(type){
    if (type == "genome"){
        'orangered4'
    }else if (type == "mitochondria"){
        'limegreen'
    }else if (type == "chloroplast"){
        'navy'
    }else if (type == '*') {
        'black'
    }else{
        print(type)
        'skyblue'
    }
}

plot_each_coverage <- function(ls){
    name <- stri_replace_all_fixed(ls$name,"/","-")
    g <- tibble(index = 1:length(ls$value),
                coverage = ls$value) %>%
        ggplot(mapping = aes(x = index,
                             y = coverage)) +
        geom_smooth(color = chrname_to_color(ls$type)) +
        cowplot::theme_cowplot() +
        labs(title = "positionwise_coverage:\nID:" %s+% ls$name %s+% "\nAssign:" %s+% ls$type)
    cowplot::ggsave(filename = "./pdf/cov_" %s+% name %s+% ".pdf")
    cowplot::ggsave(filename = "./png/cov_" %s+% name %s+% ".png")
}

positionwise_parsed %>% lapply(plot_each_coverage)
