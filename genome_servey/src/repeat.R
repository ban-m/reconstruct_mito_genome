library("tidyverse")
data <- read_tsv("./result/repeat_result.tsv",col_names = TRUE)
g <- data %>% ggplot(aes(x = frequency, y = count,color = as.factor(kmer_size)))  + geom_line()


freq_decay <- data %>% filter(frequency == 2) %>% select(count,kmer_size)
g <- freq_decay %>%
    ggplot(aes(x = kmer_size,y = count)) + geom_line()

## Seemingly, from 500?

linear <- freq_decay %>% filter(kmer_size >= 500) %>%
    (function(x)lm(count ~ kmer_size, x))


repeat_count <- function(x){
    linear$coefficients["(Intercept)"] + x * linear$coefficients["kmer_size"]
}

gg <- g + stat_function(colour = "orange",fun = repeat_count)

summary(linear)

## The coefficient is almost exactly 2, meaning only single pair of intersparsed repeat.
## As shown in "17  8389      1200", the length of that repeat is about 8389 + 1200 = 9500.

## We can eliminate this "largest" repeat just substract the expected size.
freq_decay %>% mutate(count = count - repeat_count(kmer_size))

##  8  1282          200
## Still, there are cetain amount of intersparsed repeat, causing 
