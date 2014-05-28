
#dat<-read.csv("D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset_clean.csv")

#dat<-read.csv("/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset_clean.csv")

#dat<-read.csv("/home/din02g/Google Drive/WBHC-Project/data/Train/COI_Train_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="/home/din02g/Google Drive/WBHC-Project/data/Train/COI_Train_Dataset_clean.csv")

dat<-read.table("D:/Users/Dinnage/iBOL/Insect_bold_data_2.txt",sep="\t",stringsAsFactors=FALSE,
                 header=TRUE, quote="", comment.char="")

cleandat<-dat[,c("processid","sampleid","recordID","catalognum","fieldnum","order_name",
                 "family_name","subfamily_name","genus_name","species_name","sequenceID",
                 "markercode","nucleotides")]

## get rid of data if it does not have at least a family level ID
cleandat<-cleandat[!is.na(dat$family_taxID),]

## reduce only to COI-5P barcodes

cleandat<-cleandat[cleandat$markercode=="COI-5P",]

## replace blanks with NAs
cleandat[cleandat==" "]<-NA

## remove entries with no sequence data

cleandat<-cleandat[!is.na(cleandat$nucleotides),]

## reduce to only one sequence per species
library(dplyr)

## make my own ID
cleandat$IDnum<-seq_along(cleandat$processid)

## split into smaller dataframes
cleandat.list<-split(cleandat,gl(n=ceiling(nrow(cleandat)/1000),k=1000,length=nrow(cleandat)))

path<-"D:/Users/Dinnage/Projects/WBHC-Project/data/FullData"
for (i in 1:length(cleandat.list)){
  fullpath<-paste(path,"/insect_COI_data_clean_",sprintf("%03d", i),".csv",sep="")
  write.csv(cleandat.list[[i]],fullpath,quote=FALSE,row.names=FALSE)
  print(i)
}

## make training and test set

