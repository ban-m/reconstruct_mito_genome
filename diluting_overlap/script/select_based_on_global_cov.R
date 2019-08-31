                                        # ---- import library ----

library("tidyverse")
library("stringi")
loadNamespace("cowplot")
source("~/work/generalplot.R")
                                        # ---- import arguments ----

args <- commandArgs(trailingOnly = TRUE)
## args <- c("/grid/ban-m/arabidopsis_thaliana/sequel/dilution_onthefly2/*.tsv",
##           "./result/sequel_to_reference.tsv",
##           "./result/selected.tsv")
files <- Sys.glob(args[1])
trial_num <- length(files)

                                        # --- main control flow ----
                                        # --- filtering function ----
filtering_with <- function(x,annot){
    data <- readr::read_tsv(x,col_names=FALSE) %>%
        rename(id = X1,mean=X2, sd = X3, length = X4) 
    coverage_data <- full_join(data,annot,by="id")
    cleaned_data <- coverage_data %>% filter(!is.na(mean) & sd != 0) 
    size <- cleaned_data %>% pull(length) %>% max()
    cleaned_data <- cleaned_data %>% mutate(size = length * 10 / size) %>%
        mutate(cov = sd/mean)
    cleaned_data %>% nest(-class)%>%
        mutate(data = map(data, ~ .x %>% summarize(m_mean = mean(mean), sd_mean = sd(mean),
                                                   s_mean = mean(sd), sd_sd = sd(sd),
                                                   c_mean = mean(cov), c_sd = sd(cov)))) %>%
        unnest() %>% print()
    cleaned_data %>% mutate(
                         prediction = (mean < 20) & (mean > 2.5) & (cov < 1.0) & (cov < 1 - mean / 20) ) %>%
        select(id,prediction)
}


                                        # --- parsing, opening, filtering ---
annot <- read_tsv(args[2],col_names = FALSE) %>% rename(id=X1, class=X2)
filtering <- function(x){
    filtering_with(x,annot)
}
raw_data <- files %>% map( .f = filtering) %>% reduce(.f = function(x,y){full_join(x =x,y=y,by="id")})
data <- raw_data %>% mutate(result = raw_data %>% select(-id) %>% rowSums(na.rm = TRUE)) %>% mutate( result = result/trial_num) %>% select(id,result)
selected <- data %>% filter(result > 0.5)
write_tsv(selected,args[3])

all <- annot %>% count() %>% pull(n)
mito <- annot %>% filter(class == "mitochondria") %>% count() %>% pull(n)
print("initial condition:" %s+% (1.0 * mito / all * 100))


data <- inner_join(selected,annot,by="id") 
mito_after <- data %>% filter(class == "mitochondria") %>% count() %>% pull(n)
all_after <- data %>% count() %>% pull(n)
print("After purification:" %s+% (1.0 *mito_after / all_after * 100))
print("mito:" %s+% mito %s+% "->" %s+% mito_after)

## 921736 # all
## 21444 # all mito (2%)
## 330441 # cov0
## 1104 # mito cov 0
