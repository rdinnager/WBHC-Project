## Script to load in RTD chunks for full COI sequences and concatenate them into one data file for further processing
path <- "D:/Users/Dinnage/Projects/WBHC-Project/data/FullCOI/RTDs"
numfiles <- 940
library(data.table)
library(dplyr)
library(magrittr)

getFiles <- function(i) {
  fullpath <- paste(path, "/Hexapoda_COI_segments_RTD_", sprintf("%03d",i),".csv", sep="")
  RTDs <- fread(fullpath) 
  return(RTDs)
}

RTDs <- seq_len(3) %>% data.frame %>% set_names("filenum") %>% rowwise %>% 
  do(rtds=getFiles(.$filenum)) %>% select(rtds) %>% rbindlist