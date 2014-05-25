#Script to calculate Return Time Distribution from iBOL COI Data
library(rPython)
library(data.table)
library("Biostrings")
#dat<-read.csv("D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset_clean.csv")

dat<-fread("D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset_clean.csv",stringsAsFactors=F)
seqs<-DNAStringSet(dat$nucleotides)