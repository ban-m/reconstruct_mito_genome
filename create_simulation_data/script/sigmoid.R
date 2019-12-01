library("tidyverse")


sigmoid <- function(x){
    x <- x - 50
    atan(x/5)*2/pi + 1
}

back <- ggplot(data.frame(x = c(0,70)),aes(x))

g <- back + stat_function(fun = sigmoid)
