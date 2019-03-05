library("tidyverse")
source("~/work/generalplot.R")
library("stringi")


files <- c("./result/ont_32mer_freq.tsv","./result/sequel_32mer_freq.tsv")

main <- function(file){
### ===== Open File =====
    data <- read_delim(file,col_names=FALSE,delim=' ',col_types = "ii")
    outname <- basename(tools::file_path_sans_ext(file))
    g <- data %>% ggplot(mapping = aes(x = X1, y=X2)) +
        geom_point()  + 
        ylim(c(0,1000)) # +
    ## xlim(c(0,2500))
    generalplot(g,outname)
}

plot_kmer <- function(kmer){
### ---- specify file ----
    file <- "./result/sequel_" %s+% kmer %s+% "mer_freq.tsv"
    main(file)
}


mer <- c(8,9,10,12)
mer %>% sapply(plot_kmer)

files %>% sapply(main)
