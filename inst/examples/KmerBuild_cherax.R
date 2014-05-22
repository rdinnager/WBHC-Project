setwd("/data/noflush/din02g")

library("Biostrings")
library("plyr")
library("doSNOW")

source("RTD.R")

dat<-read.csv("COI_Play_Dataset.csv",stringsAsFactors=F)

dna<-dat$nucleotides
names(dna)<-c(1:length(dna))

ibolSeq<-DNAStringSet(dna)

#Generate Return Time Distribution for all COI sequences
# Generate Kmer names
maxK<-7 #maximum K-value
alph<-alphabet(ibolSeq,baseOnly=T) #DNA alphabet
kmernames<-do.call(c,lapply(c(1:maxK),function(x) mkAllStrings(alph,x)))
names(kmernames)<-c(1:length(kmernames))

#Get RTDS in parallel
#setDefaultClusterOptions(outfile="")
#setDefaultClusterOptions(outfile="/dev/null")
cl <- makeCluster(16, "MPI")
clusterExport(cl, c("get_RTDs","get_RTD","get_RTs"))
registerDoSNOW(cl)

runRTD<-function(x,data) {
  result<-get_RTDs(x,data)
  print(as.numeric(names(x)))
  return(result)
}
RTDs<-ldply(kmernames,get_RTDs,data=ibolSeq,.parallel=T,.paropts=list(.export=c("ibolSeq"),.packages="Biostrings"))
stopCluster(cl)

RTDs<-t(as.matrix(RTDs[1:(2*length(kmernames)),2:(ncol(RTDs))]))
RTDs[is.na(RTDs)]<-0
knames<-as.vector(rbind(paste(kmernames,"_Mean",sep=""),paste(kmernames,"_Var",sep="")))
colnames(RTDs)<-knames

RTDnames<-dat[,c("processid","order_name","family_name","subfamily_name","genus_name","species_name")]

write.csv(RTDs,file="Play_RTDs.csv",row.names=F)
write.csv(RTDnames,file="Play_Groups.csv",row.names=F)