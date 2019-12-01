library("tidyverse")
loadNamespace("cowplot")
## ===== Setting for usual use ======
generalplot <- function(g,name){
    cowplot::ggsave(filename = paste0("./pdf/",name,".pdf"),
                    plot = g + cowplot::theme_cowplot())
    cowplot::ggsave(filename = paste0("./png/",name,".png"),
                    plot = g + cowplot::theme_cowplot())
}

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
outputname <- args[2]

dataset <- read_tsv(filename)

g <- dataset %>% ggplot() +
    geom_point(mapping = aes(x = TP_EM, y = TN_EM, color = Skew, size = Coverage),alpha = 0.4) +
    scale_colour_gradient2(mid = "black",low = "red", high = "blue", midpoint = 0.5)
generalplot(g,paste0(outputname,"_EM_point"))


g <- dataset %>% ggplot() +
    geom_point(mapping = aes(x = TP_Naive, y = TN_Naive, color = Skew, size = Coverage,alpha = 0.4)) +
    scale_colour_gradient2(mid = "black",low = "red", high = "blue", midpoint = 0.5)
generalplot(g,paste0(outputname,"_Naive_point"))

g <- dataset %>% mutate(TP_Improve = TP_EM-TP_Naive,
                        TN_Improve = TN_EM-TN_Naive) %>%
    ggplot() +
    geom_point(mapping = aes(x = TP_Improve, y = TN_Improve, color = Skew, size = Coverage,alpha = 0.4)) +
    scale_colour_gradient2(mid = "black",low = "red", high = "blue", midpoint = 0.5)
generalplot(g,paste0(outputname,"_Diff_point"))


g <- dataset %>% mutate(TP_Improve = TP_EM-TP_Naive,
                        TN_Improve = TN_EM-TN_Naive) %>%
    select(TP_Improve, TN_Improve, Skew, Coverage) %>%
    nest(-Skew, -Coverage) %>%
    mutate(data = map(data, function(x) x %>% summarize(TP_Improve_mean = mean(TP_Improve), TN_Improve_mean = mean(TN_Improve))))%>%
    unnest() %>% 
    ggplot() +
    geom_point(mapping = aes(x = TP_Improve_mean,
                             y = TN_Improve_mean,
                             color = Skew, size = Coverage,alpha = 0.4)) +
    scale_colour_gradient2(mid = "black",low = "red", high = "blue", midpoint = 0.5)
generalplot(g,paste0(outputname,"_Diff_point_mean"))


fscores <- dataset %>%
    mutate(EM_f = 2 * TP_EM * TN_EM / (TP_EM + TN_EM),
           Naive_f = 2 * TP_Naive * TN_Naive / (TP_Naive + TN_Naive)) %>%
    mutate(f_diff = EM_f - Naive_f) %>% 
    select(Coverage, Skew, EM_f, Naive_f, f_diff)


g <- fscores %>%
    ggplot() +
    geom_violin(aes(x = factor(Skew), y = f_diff)) + facet_wrap(.~Coverage)
generalplot(g,paste0(outputname,"_fscores_diff_violin"))


g <- fscores %>%
    nest(-Coverage, -Skew) %>% 
    mutate(data = map(data, function(x) x %>% summarize(f_diff= mean(f_diff,na.rm = TRUE ))))%>%
    unnest() %>% 
    ggplot() +
    geom_line(mapping = aes(x = Skew, y = f_diff))

