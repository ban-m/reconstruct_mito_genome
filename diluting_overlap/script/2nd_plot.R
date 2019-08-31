                                        # ---- import library ----

library("tidyverse")
library("stringi")
loadNamespace("cowplot")
source("~/work/generalplot.R")
                                        # ---- import arguments ----

args <- commandArgs(trailingOnly = TRUE)
## args <- c("./result/remap.tsv",
##           "./result/sequel_to_reference.tsv",
##           "./result/selected_2nd.tsv")
data <- read_tsv(args[1],col_names=FALSE) %>%
    rename(mean=X2, sd = X3, length = X4, id = X1)
annot <- read_tsv(args[2],col_names = FALSE) %>% rename(id=X1, class=X2)
coverage_data <- full_join(data,annot) 
cleaned_data <- coverage_data %>% filter(!is.na(mean))
size <- cleaned_data %>% pull(length) %>% max()
cleaned_data <- cleaned_data %>% mutate(size = length * 10 / size) %>%
    mutate(cov = sd/mean)
summary <- cleaned_data %>% nest(-class)%>%
    mutate(data = map(data, ~ .x %>% summarize(m_mean = mean(mean), sd_mean = sd(mean),
                                               s_mean = mean(sd), sd_sd = sd(sd),
                                              c_mean = mean(cov), c_sd = sd(cov)))) %>%
    unnest()
write_tsv(x = summary, path = "./result/2dn_summary.tsv")


g <- cleaned_data %>% sample_n(10000) %>%
    ggplot(mapping = aes(x = mean, y = sd, size = size,color = class)) + geom_point(alpha=0.4)
generalplot(g,"2nd_all")
cleaned_data %>% filter(mean > 200) %>% write_tsv(path = args[3])
