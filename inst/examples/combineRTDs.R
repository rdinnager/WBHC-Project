## Script to load in RTD chunks for full COI sequences and concatenate them into one data file for further processing
path <- "D:/Users/Dinnage/Projects/WBHC-Project/data/FullCOI/RTDs"
numfiles <- 940
library(data.table)
library(dplyr)
library(magrittr)

types <- read.csv(paste(path, "/Hexapoda_COI_segments_RTD_", sprintf("%03d",1),".csv", sep=""), nrows=2, 
                          stringsAsFactors=FALSE, header=TRUE)[2,]
types[7:10926]<-types[7:10926]+0.0123456789
types<-as.list(types)
getFiles <- function(i) {
  fullpath <- paste(path, "/Hexapoda_COI_segments_RTD_", sprintf("%03d",i),".csv", sep="")
  RTDs <- scan(fullpath, what = types, sep=",", skip = 1) 
  attr(RTDs, "row.names") <- .set_row_names(length(RTDs[[1]]))
  class(RTDs) <- "data.frame"
  return(RTDs)
}

RTDs <- seq_len(numfiles) %>% data.frame %>% set_names("filenum") %>% rowwise %>% 
  do(getFiles(.$filenum))