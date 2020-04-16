library("tidyverse")
loadNamespace("cowplot")
library(Rtsne)
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi=350,width = 178,height = 86,units="mm")
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi = 350,width = 178,height = 86,units="mm")
    ## cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
    ##                 plot = g + cowplot::theme_cowplot())
    ## cowplot::ggsave(filename = paste0("./png/",name,".png"),
    ##                 plot = g + cowplot::theme_cowplot())
}

### ---- Two cluster -----
dataset <- read_tsv("./result/lks_two.tsv",col_names=FALSE) %>%
    rename(State = X2, ID=X3, ReadID =X4, Cluster=X5, Prob1= X6, Prob2=X7) %>%
    mutate(Cluster = factor(Cluster)) %>%
    mutate(State = ifelse(State == "B", "BEFORE", "AFTER"))

g <- dataset %>%
    ggplot() + geom_point(mapping = aes(x = Prob1, y = Prob2, color = Cluster)) +
    facet_grid(State~.) +
    labs(x = "P(Read|Model1)", y = "P(Read|Model0)")
generalplot(g,"two_clusters")



### ---- Six cluster -----
dataset <- read_tsv("./result/lks_six.tsv",col_names=FALSE) %>%
    rename(State = X2, ID=X3, ReadID =X4, Cluster=X5) %>%
    rename(LK1= X6, LK2=X7, LK3=X8, LK4=X9,LK5 = X10,LK6 = X11) %>%
    mutate(State = ifelse(State == "B", "BEFORE", "AFTER")) %>% 
    select(-ID, -ReadID, -X1)

initial_state <- dataset %>% filter(State == "BEFORE") %>% select(-State, -Cluster)
answer <- dataset %>% filter(State == "BEFORE") %>% pull(Cluster)
result <- Rtsne(X=as.matrix(initial_state),pca=FALSE,theta=0,perplexity = 100,  normalize=TRUE,max_iter = 10000)
is_tsne <- tibble(tSNE1= result$Y[,1],tSNE2=result$Y[,2], answer = answer, State = "BEFORE")
final_state <- dataset %>% filter(State == "AFTER") %>% select(-State, -Cluster)
answer <- dataset %>% filter(State == "AFTER") %>% pull(Cluster)
result <- Rtsne(X=as.matrix(final_state),pca=FALSE,theta=0, perplexity = 100,normalize=TRUE, max_iter = 10000)
fs_tsne <- tibble(tSNE1= result$Y[,1],tSNE2=result$Y[,2], answer = answer, State= "AFTER")
tsne_data <- bind_rows(is_tsne, fs_tsne) 
g <- tsne_data %>% mutate(Cluster = factor(answer)) %>% 
    ggplot() +
    geom_point(aes(x = tSNE1,y = tSNE2, color = Cluster)) +
    facet_grid(State~.)

generalplot(g,"six_clusters")



