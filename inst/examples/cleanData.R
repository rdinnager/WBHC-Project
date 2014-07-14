
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

## remove entries with no sequence data or very small sequences

cleandat<-cleandat[!is.na(cleandat$nucleotides),]
cleandat<-cleandat[nchar(cleandat$nucleotides)>50,]

## Strip out gap characters
cleandat$nucleotides <- gsub("-", "", cleandat$nucleotides, fixed = TRUE)

## reduce to only one sequence per species
library(dplyr)
## data with species names
cleandat.spec<-cleandat[!is.na(cleandat$species_name),]
#cleandat.new<-group_by(cleandat.spec,species_name) %>% do(sample_n(.,1))
## take the sequence with the most number of bases for each duplicated species
cleandat.new<-group_by(cleandat.spec,species_name) %>% do(.[which.max(nchar(.$nucleotides)),])
cleandat<-cleandat[is.na(cleandat$species_name),]
cleandat<-rbind(cleandat,cleandat.new)

## make my own ID
cleandat$IDnum<-seq_along(cleandat$processid)

## randomize row order
cleandat<-cleandat[sample.int(nrow(cleandat)),]

#collect some data for later use
sumdat<-group_by(cleandat,order_name,family_name,genus_name) %>% summarise(count=n())
IDnums<-group_by(cleandat,order_name,family_name,genus_name) %>% do(IDnums=.$IDnum)
fullsumdat<-left_join(sumdat,IDnums)

## save data
path<-"D:/Users/Dinnage/Projects/WBHC-Project/data/FullData"
#save(fullsumdat,file=paste(path,"/insect_COI_data_summary.Rdata",sep=""))
saveRDS(fullsumdat,file=paste(path,"/insect_COI_data_summary_July14_2014.rds",sep=""))
saveRDS(cleandat,file=paste(path,"/insect_COI_data_full_July14_2014.rds",sep=""))
#write.csv(sumdat,file=paste(path,"/insect_COI_data_counts.csv",sep=""), row.names=FALSE)
write.csv(cleandat,file=paste(path,"/insect_COI_data_full_July14_2014.csv",sep=""), row.names=FALSE)
#dput(fullsumdat,file=paste(path,"/insect_COI_data_sumtest.txt",sep=""))

## write sequences to fasta file
library(Biostrings)
seqs <- DNAStringSet(cleandat$nucleotides)
names(seqs) <- cleandat$IDnum
testset <- seqs[1:10]
writeXStringSet(seqs, paste(path,"/insect_COI_data_sequences_July14_2014.fasta",sep=""))
writeXStringSet(testset, paste(path,"/insect_COI_data_forTesting_July14_2014.fasta",sep=""))

## split into smaller dataframes of 500 rows each
cleandat.list<-split(cleandat,gl(n=ceiling(nrow(cleandat)/500),k=500,length=nrow(cleandat)))

for (i in 1:length(cleandat.list)){
  fullpath<-paste(path,"/insect_COI_data_clean_",sprintf("%03d", i),".csv",sep="")
  write.csv(cleandat.list[[i]],fullpath, row.names=FALSE)
  print(i)
}

## make training and test set

