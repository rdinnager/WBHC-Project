#Script to calculate Return Time Distribution from iBOL COI Data
library(rPython)
library(data.table)
library("Biostrings")
#dat<-read.csv("D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset_clean.csv")
#dat<-read.csv("/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset_clean.csv")

#dat<-fread("D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset_clean.csv",stringsAsFactors=F)
dat<-fread("/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset_clean.csv",stringsAsFactors=F)
## get rid of empty sequences
dat<-dat[dat$nucleotides!="",]
## convert nucleotides to Biostrings sequence object
seqs<-DNAStringSet(dat$nucleotides)

test<-seqs
names(test)<-seq_along(test)
writeXStringSet(test,"/home/din02g/Google Drive/WBHC-Project/data/testing.fa")

tester<-getRTDs("/home/din02g/Google Drive/WBHC-Project/data/testing.fa",k=1)
system.time(
testall<-getAll_RTDs("/home/din02g/Google Drive/WBHC-Project/data/testing.fa",k=7)
)