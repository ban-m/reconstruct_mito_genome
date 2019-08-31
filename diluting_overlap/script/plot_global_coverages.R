                                        # ---- import library ----

library("tidyverse")
library("stringi")
loadNamespace("cowplot")
source("~/work/generalplot.R")
                                        # ---- import arguments ----

args <- commandArgs(trailingOnly = TRUE)
## args <- c("/grid/ban-m/arabidopsis_thaliana/sequel/dilution_onthefly/dilute.tsv",
##           "./result/sequel_to_reference.tsv"
##           "./result/selected.tsv")
data <- readr::read_tsv(args[1],col_names=FALSE) %>%
    rename(mean=X2, sd = X3, length = X4, id = X1)
annot <- read_tsv(args[2],col_names = FALSE) %>% rename(id=X1, class=X2)
coverage_data <- full_join(data,annot) 
cleaned_data <- coverage_data %>% filter(!is.na(mean) & sd != 0) 
size <- cleaned_data %>% pull(length) %>% max()
cleaned_data <- cleaned_data %>% mutate(size = length * 10 / size) %>%
    mutate(cov = sd/mean)


summary <- cleaned_data %>% nest(-class)%>%
    mutate(data = map(data, ~ .x %>% summarize(m_mean = mean(mean), sd_mean = sd(mean),
                                               s_mean = mean(sd), sd_sd = sd(sd),
                                               c_mean = mean(cov), c_sd = sd(cov)))) %>%
    unnest()

write_tsv(x = summary, path = "./result/summary.tsv")


all <- annot %>% count() %>% pull(n)
mito <- annot %>% filter(class == "mitochondria") %>% count() %>% pull(n)
print("initial condition:" %s+% (1.0 * mito / all * 100))
                                        # --- main control flow ----
g <- cleaned_data %>% sample_n(10000) %>%
    ggplot(mapping = aes(x = mean, y = sd, size = size,color = class)) + geom_point(alpha=0.4)
generalplot(g,"all")
g <- cleaned_data %>% sample_n(10000) %>%
    ggplot(mapping = aes(x = mean, y = cov, size = size,color = class))  + geom_point()
generalplot(g,"all_cov")

g <- cleaned_data %>% filter(mean < 20 & 2 < mean) %>% 
    ggplot(mapping = aes(x = mean, y = sd, size = size,color = class))  + geom_point()
generalplot(g,"all_mean_20")
g <- cleaned_data %>% filter(mean < 20 & 2 < mean) %>% 
    ggplot(mapping = aes(x = mean, y = cov, size = size,color = class))  + geom_point()
generalplot(g,"all_cov_mean_20")

mitohex <-  cleaned_data %>% filter(class == "mitochondria") %>%
    filter(mean < 20) %>% filter(mean > 5) %>% 
    ggplot(mapping = aes(x = mean, y = sd)) + geom_hex() 
generalplot(mitohex,"mito_hex")

hex <- cleaned_data %>% sample_n(10000) %>%
    filter(mean < 20) %>% filter(mean > 5) %>% 
    ggplot(mapping = aes(x = mean, y = sd)) + geom_hex()
generalplot(hex,"hex_all")

data <- cleaned_data %>% filter(mean < 20) %>% filter(mean > 2.5) %>% filter(cov < 1.0) %>% filter(cov < 1 - mean / 20)
g <- data %>% ggplot(mapping = aes(x = mean, y = cov, color = class)) + geom_point(alpha = 0.4)
generalplot(g,"result")
g <- data %>% filter(class == "mitochondria") %>% ggplot(mapping = aes(x = mean, y = cov, color = class)) + geom_point(alpha = 0.4)
generalplot(g,"result_mito")
write_tsv(x = data,path = args[3])

all_after <- data %>% count() %>% pull(n)
mito_after <- data %>% filter(class == "mitochondria") %>% count() %>% pull(n)
print("After purification:" %s+% (1.0 *mito_after / all_after * 100))
print("mito:" %s+% mito %s+% "->" %s+% mito_after)

## 921736 # all
## 21444 # all mito (2%)
## 330441 # cov0
## 1104 # mito cov 0
