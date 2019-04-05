library("tidyverse")
library("stringi")
source("~/work/generalplot.R")
args <- commandArgs(trailingOnly = TRUE)

#### ===== Load data ======

args <- c("./result/sequel_positionwise_start.pick.csv",
          "./result/sequel_positionwise_stop.pick.csv",
          "./result/sequel_minimap2_assignment.tsv")

assignment <- read_tsv(args[3],col_name = FALSE, col_types = "cc") %>%
    rename(id = X1, type = X2)
assign_table <- assignment %>% mutate(chrtype = ifelse(stri_detect_charclass(type,"[0-9]"),"genome",type)) %>%
    pull(chrtype)
names(assign_table) <- assignment$id


to_vectors <- function(input){
    strsplit(x = input, split = ',') %>% lapply(function(ls){
        name <- ls[1]
        datasize <- as.numeric(ls[2])
        position <- as.numeric(ls[3:(2+datasize)])
        value <- as.numeric(ls[(3+datasize):(2+datasize*2)])
        list(name = name,
             position = position,
             type = assign_table[name],
             value = value,
             datasize = datasize
             )
    })
}


open_position_file <- function(file){
    to_vectors(readLines(file))
}

starts <- open_position_file(args[1])
stops <- open_position_file(args[2])

marge_two <- function(start_data,stop_data){
    print(start_data$name %s==% stop_data$name)
    print(start_data$datasize)
    data <- tibble(position = c(start_data$position,stop_data$position),
                   count = c(start_data$value, stop_data$value),
                   type = c(rep("start",start_data$datasize),rep("end",stop_data$datasize)))
    list(data = data,
         name = start_data$name)
}

read_informations <- map2(starts, stops,marge_two) 

plot_read <- function(ls){
    name <- stri_replace_all_fixed(ls$name,"/","-")
    g <- ls$data %>% ggplot(mapping = aes(x = position , y = count,  color = type, fill = type)) +
        geom_point() +
        labs(title = "ReadID:" %s+% ls$name %s+% "\nAssing:" %s+% ls$type,
             x = "position of read" , y = "read count")
    generalplot(g=g,name="stst_" %s+% name)
}

result <- map(read_informations,plot_read)
