library("tidyverse")
canu <- read_tsv("./result/canu_diff.txt")
canu_summary <- canu %>%
    summarize(both = sum(IsInFirst & IsInSecond),
              only_bootstrap = sum(IsInFirst) - both,
              only_mapback = sum(IsInSecond) - both)
write_tsv(canu_summary, "./result/canu_diff.summary")

flye <- read_tsv("./result/flye_diff.txt")
flye_summary <- flye %>% 
    summarize(both = sum(IsInFirst & IsInSecond),
              only_bootstrap = sum(IsInFirst) - both,
              only_mapback = sum(IsInSecond) - both)
write_tsv(flye_summary, "./result/flye_diff.summary")


wtdbg <- read_tsv("./result/wtdbg_diff.txt")
wtdbg_summary <- wtdbg %>% 
    summarize(both = sum(IsInFirst & IsInSecond),
              only_bootstrap = sum(IsInFirst) - both,
              only_mapback = sum(IsInSecond) - both)
write_tsv(wtdbg_summary, "./result/wtdbg_diff.summary")
