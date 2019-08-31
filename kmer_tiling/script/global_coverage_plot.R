library("tidyverse")
library("stringi")
loadNamespace("cowplot")
source("~/work/generalplot.R")
args <- commandArgs(trailingOnly = TRUE)

### ------ Load data -----

## args <- c("./result/tiling_20.tsv",
##           "./result/sequel_minimap2_assignment.tsv")

data <- read_tsv(args[1])
max_len <- data %>% pull(Len) %>% max();
annot <- read_tsv(args[2],col_names = FALSE,col_types = c("cc")) %>%
    rename(ID = X1, type = X2)
data <- full_join(data,annot) %>% filter(!is.na(Mean))  %>%
    mutate(chrtype = ifelse(stri_detect_charclass(type,"[0-9]"),"genome",type)) %>%
    mutate(size = Len * 10 / max_len )
filename <- stri_replace_all_fixed(basename(args[1]),".tsv","")

makeFileName <- function(filename,adress){
    c("./pdf/" %s+% filename %s+% adress %s+% ".pdf",
      "./png/" %s+% filename %s+% adress %s+% ".png")
}

g <- data %>% sample_n(10000) %>% ggplot(aes(x = SD, y = Mean, size = size, color = chrtype)) +
    geom_point(alpha = 0.4)
makeFileName(filename,"_all") %>% sapply(FUN=function(filename)cowplot::ggsave(filename = filename,plot=g))
g <- data %>% filter(type == "mitochondria") %>% ggplot(aes(x = SD, y = Mean, size = size, color = chrtype)) +
    geom_point()
makeFileName(filename,"_mito") %>% sapply(FUN=function(filename)cowplot::ggsave(filename = filename,plot=g))
