library(tidyverse)
loadNamespace("cowplot")
library(Rtsne)
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi=600,width = 178,height = 86,units="mm")
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot(font_size=12),
                    dpi = 600,width = 178,height = 86,units="mm")
}


## ===== Plot posterior probability plots ================
### ---- Two cluster -----
dataset <- read_tsv("./result/posterior_probability/two_cluster_posterior.tsv",
                    col_names=FALSE) %>%
    rename(State = X1, ID=X2, ReadID =X3, Cluster=X4, Prob1= X5, Prob2=X6) %>%
    mutate(Cluster = factor(Cluster)) %>%
    mutate(State = ifelse(State == "B", "BEFORE", "AFTER"))

g <- dataset %>%
    ggplot() + geom_point(mapping = aes(x = Prob1, y = Prob2, color = Cluster)) +
    facet_grid(State~.) +
    labs(x = "P(Read|Model 0)", y = "P(Read|Model 1)")
generalplot(g,"two_clusters")


### ---- Six cluster -----
dataset <- read_tsv("./result/posterior_probability/six_cluster_posterior.tsv",
                    col_names=FALSE) %>%
    rename(State = X1, ID=X2, ReadID =X3, Cluster=X4) %>%
    rename(LK1= X5, LK2=X6, LK3=X7, LK4=X8,LK5 = X9,LK6 = X10) %>%
    mutate(State = ifelse(State == "B", "BEFORE", "AFTER")) %>% 
    select(-ID, -ReadID)

initial_state <- dataset %>% filter(State == "BEFORE") %>% select(-State, -Cluster)
answer <- dataset %>% filter(State == "BEFORE") %>% pull(Cluster)
result <- Rtsne(X=as.matrix(initial_state),pca=FALSE,theta=0,perplexity = 100,  normalize=TRUE,max_iter = 10000)
is_tsne <- tibble(tSNE1= result$Y[,1],tSNE2=result$Y[,2], answer = answer, State = "BEFORE")
final_state <- dataset %>% filter(State == "AFTER") %>% select(-State, -Cluster)
answer <- dataset %>% filter(State == "AFTER") %>% pull(Cluster)
result <- Rtsne(X=as.matrix(final_state),pca=FALSE,theta=0, perplexity=100,normalize=TRUE, max_iter = 10000)
fs_tsne <- tibble(tSNE1= result$Y[,1],tSNE2=result$Y[,2], answer = answer, State= "AFTER")
tsne_data <- bind_rows(is_tsne, fs_tsne) 

g <- tsne_data %>% mutate(Cluster = factor(answer)) %>% 
    ggplot() +
    geom_point(aes(x = tSNE1,y = tSNE2, color = Cluster)) +
    facet_grid(State~.)

generalplot(g,"six_clusters")



## ========= Plot the result of the benchmarks ============

dataset <- read_tsv("./result/benchmark/benchmark.tsv", col_names=FALSE)
outname <- "benchmark"
accs <- dataset %>% mutate(acc1 = X4 + X7, acc2 = X5 + X6, acc = pmax(acc1, acc2)) %>%
    mutate(acc = acc / X9) %>% select(X2,X3,acc1,acc2, acc, X8,X9,X10, X11)

g <- accs %>%
    filter(X10 <= 4) %>%
    ggplot() + geom_boxplot(mapping = aes(x = factor(X9), y = acc, color = factor(X10))) +
    labs(x = "Total Coverage", y = "Accuracy", color = "Number of Variants\nout of 9K bp")

generalplot(g,paste0(outname,"_boxplot"))

g <- accs %>% filter(X9 %in% c(80, 110, 140)) %>%
    ggplot() + geom_boxplot(mapping = aes(x = factor(X9), y = acc, color = factor(X11))) +
    labs(x = "Total Coverage", y = "Accuracy", color = "Number of Variants\nout of 9K bp") +
    coord_flip() + facet_grid(.~factor(X10))

generalplot(g,paste0(outname,"_boxplot_flip"))
