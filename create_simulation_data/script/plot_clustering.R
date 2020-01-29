library("tidyverse")
loadNamespace("cowplot")
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
dataset <- read_tsv("./result/lks.tsv",col_names=FALSE)
g <- dataset %>% nest(-X4) %>% head(n=1) %>% unnest() %>%
    mutate(answer = factor(as.integer(X3 < 50)))%>%
    rename(LK1 = X1, LK2 = X2) %>% 
    ggplot() + geom_point(aes(x = LK1, y = LK2, color = answer))
generalplot(g,"2_initial")

g <- dataset %>% nest(-X4) %>% tail(n=1) %>%
    unnest() %>%
    mutate(answer = factor(as.integer(X3 < 50)))%>%
    rename(LK1 = X1, LK2 = X2) %>% 
    ggplot() + geom_point(aes(x = LK1, y = LK2, color = answer))
generalplot(g,"2_final")

diffs <- dataset %>% nest(-X4) %>%
    mutate(data = map(data, ~ summarize(.,diff = sum(abs(X1 - X2))))) %>%
    unnest()
g <- diffs %>% ggplot() + geom_point(aes(x = X4, y = diff))    


### ---- Two cluster -----
dataset <- read_tsv("./result/lks_skew.tsv",col_names=FALSE)
g <- dataset %>% nest(-X4) %>% head(n=1) %>% unnest() %>%
    mutate(answer = factor(as.integer(X3 < 30))) %>%
    rename(LK1 = X1, LK2 =X2) %>% 
    ggplot() + geom_point(aes(x = LK1, y = LK2, color = answer))
generalplot(g,"2_skew_initial")

g <- dataset %>% nest(-X4) %>% tail(n=1) %>%
    unnest() %>%
    mutate(answer = factor(as.integer(X3 < 30))) %>%
    rename(LK1 = X1, LK2 =X2) %>% 
    ggplot() + geom_point(aes(x = LK1, y = LK2, color = answer))
generalplot(g,"2_skew_final")

diffs <- dataset %>% nest(-X4) %>%
    mutate(data = map(data, ~ summarize(.,diff = sum(abs(X1 - X2))))) %>%
    unnest()
g <- diffs %>% ggplot() + geom_point(aes(x = X4, y = diff))    


### ---- Six cluster -----
dataset <- read_tsv("./result/lks_multi.tsv",col_names=FALSE)
answer <- c(rep(1:6,each =2), rep(1:6,each=49))
initial_state <- dataset %>% nest(-X8) %>% head(n=1) %>% unnest() %>%
    select(-X8, -X7)
result <- Rtsne(X=as.matrix(initial_state),pca=FALSE,theta=0, normalize=FALSE,eta = 100)
is_tsne <- tibble(tSNE1= result$Y[,1],tSNE2=result$Y[,2], answer = factor(answer))
g <- is_tsne %>%
    ggplot() +
    geom_point(aes(x = tSNE1,y = tSNE2, color = factor(answer)))
generalplot(g,"6_initial")

final_state <- dataset %>% nest(-X8) %>% tail(n=1) %>% unnest() %>%
    select(-X8,-X7)
result <- Rtsne(X=as.matrix(final_state),pca=FALSE,theta=0, normalize=FALSE,eta = 100)
fs_tsne <- tibble(tSNE1= result$Y[,1],tSNE2=result$Y[,2], answer = factor(answer))
g <- fs_tsne %>%
    ggplot() +
    geom_point(aes(x = tSNE1,y = tSNE2, color = factor(answer)))
generalplot(g,"6_final")



