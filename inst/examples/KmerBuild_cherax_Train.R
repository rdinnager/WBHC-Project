setwd("/data/noflush/din02g")
#setwd("/home/din02g/Kmer Hierarchical Classification")

library("Biostrings")
library("plyr")
library("doSNOW")

source("RTD.R")

dat<-read.csv("COI_Train_Dataset.csv",stringsAsFactors=F)

dna<-dat$nucleotides
names(dna)<-c(1:length(dna))

ibolSeq<-DNAStringSet(dna)

splits<-ceiling(length(ibolSeq)/5000)

gs<-gl(splits,5000,length(ibolSeq))
SeqList<-split(ibolSeq,gs)

datList<-split(dat,gs)


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

for (i in 1:splits){
  Seqs<-SeqList[[i]]
  RTDs<-ldply(kmernames,get_RTDs,data=Seqs,.parallel=T,.paropts=list(.export=c("Seqs"),.packages="Biostrings"))
  

  RTDs<-t(as.matrix(RTDs[1:(2*length(kmernames)),2:(ncol(RTDs))]))
  RTDs[is.na(RTDs)]<-0
  knames<-as.vector(rbind(paste(kmernames,"_Mean",sep=""),paste(kmernames,"_Var",sep="")))
  colnames(RTDs)<-knames

  RTDnames<-datList[[i]][,c("processid","order_name","family_name","subfamily_name","genus_name","species_name")]

  write.csv(RTDs,file=paste("Train_RTDs_",i,".csv",sep=""),row.names=F)
  write.csv(RTDnames,file=paste("Train_Groups_",i,".csv",sep=""),row.names=F)
  print(i)
}
stopCluster(cl)
