library("tidyverse")
library("stringi")
source("~/work/generalplot.R")
args <- commandArgs(trailingOnly = TRUE)

### ----- Load Data -----
args <- c("./result/sequel_medmid_coverage.tsv",
          "./result/sequel_minimap2_assignment.tsv")

data <- read_tsv(args[1],col_names = FALSE) %>%
    rename(id = X1, median = X2, mode = X3, length = X4)
annot <- read_tsv(args[2],col_names = FALSE,col_types = c("cc")) %>%
    rename(id = X1, type = X2)
data <- full_join(data,annot) %>% filter(!is.na(median))  %>%
    mutate(chrtype = ifelse(stri_detect_charclass(type,"[0-9]"),"genome",type)) 


gm <- data %>% filter(chrtype != "chloroplast") %>%
    ggplot(mapping = aes(x=median, y = mode,fill = chrtype, color = chrtype))
gh <-  data %>% ggplot(mapping = aes(x=median,fill = chrtype, color = chrtype))
gl <-  data %>% ggplot(mapping = aes(x=mode,fill = chrtype, color = chrtype))

g <- gh + 
    geom_histogram(bins = 60, alpha = 0.3, position = "identity") 
generalplot(g = g, name = "sequel_cov_median_overview")

g <- gh +
    geom_histogram(bins = 60, alpha = 0.3, position = "identity") +
    xlim(c(70,1000)) 
generalplot(g = g, name = "sequel_cov_median_zoom1")

g <- gh +
    geom_histogram(bins = 60, alpha = 0.3, position = "identity") +
    xlim(c(100,750)) 
generalplot(g = g, name = "sequel_cov_median_zoom2")


g <- gl +
    geom_histogram(bins = 60, alpha = 0.3, position = "identity")
generalplot(g = g, name = "sequel_cov_mode_overview")


g <- gl +
    geom_histogram(bins = 60, alpha = 0.3, position = "identity") +
    xlim(c(70,1000))
generalplot(g = g, name = "sequel_cov_mode_zoom1")

g <- gl +
    geom_histogram(bins = 60, alpha = 0.3, position = "identity") +
    xlim(c(100,750)) 
generalplot(g = g, name = "sequel_cov_mode_zoom2")


g <-  data %>% filter(chrtype == "genome") %>% ggplot(mapping = aes(x=mode,fill = chrtype, color = chrtype)) +
    geom_histogram(bins = 60, alpha = 0.3, position = "identity") +
    xlim(c(70,1000))
generalplot(g = g, name = "sequel_cov_mode_onlygenome")


g <- gm + geom_point(alpha = 0.3) 
generalplot(g = g, name = "sequel_cov_mode_median_overview")

g <- gm + geom_point(alpha = 0.3) +
    xlim(c(50,1000)) + ylim(c(50,1000))
generalplot(g = g, name = "sequel_cov_mode_median_zoom")
