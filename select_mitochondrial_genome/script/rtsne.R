library("tidyverse")
library("stringi")
loadNamespace("cowplot")
library("Rtsne")

args <- c("./result/sequel_positionwise_coverage.norm.tsv",
          "./result/sequel_minimap2_assignment.tsv")

positionwise_coverage <- read_csv(args[1],col_names =FALSE)
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

positionwise_coverage <- positionwise_coverage %>% mutate(chrtype = assign_table[X1])
subsample <- positionwise_coverage %>% sample_frac(0.05)
tsne <- Rtsne(X=subsample %>% select(-X1,-chrtype))

data <- tibble(X1 = tsne$Y[,1],X2=tsne$Y[,2],chrtype = subsample$chrtype)
g <- data %>% ggplot(mapping = aes(x = X1, y = X2, color = chrtype)) +
    geom_point(alpha = 0.3) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/rtsne.pdf",plot = g)
cowplot::ggsave(filename = "./png/rtsne.png",plot = g)

pca <- prcomp(x=subsample %>% select(-X1,-chrtype))
data <- tibble(x = pca$x[,1],y=pca$x[,2],chrtype = subsample$chrtype)
g <- data %>% ggplot(aes(x = x ,y = y , color = chrtype)) + geom_point()  +
    geom_point(alpha = 0.3) +
    cowplot::theme_cowplot()
cowplot::ggsave(filename = "./pdf/pcs.pdf",plot = g)
cowplot::ggsave(filename = "./png/pcs.png",plot = g)

