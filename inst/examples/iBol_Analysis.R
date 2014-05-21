setwd("D:/Users/Dinnage/iBol")
#ibol<-read.delim("Insect_bold_data.txt")
#ibol<-ibol[!is.na(ibol$species_taxID),]
rownames(ibol)<-ibol$recordID

#ibol$nucleotides<-sub(" ","-",ibol$nucleotides,fixed=T)
#write.csv(ibol,file="ibol_insect_species.csv")

setwd("D:/Users/Dinnage")
ibol<-read.csv("ibol_insect_species.csv")

G<-factor(ibol$genus_name)
genera<-as.character(levels(G))
write.csv(genera,file="iBolgenera.csv")
genera<-read.csv("iBolgenera.csv")
genera<-genera[-1]
library(taxize)
library(plyr)

rnames<-gnr_resolve(genera,resolve_once=T,with_context=T,http="post")
rnames<-tnrs(genera,getpost="POST",source_="NCBI",splitby=100)

ibol$nucleotides<-gsub("-","",ibol$nucleotides,fixed=T)
counts<-table(ibol$order_name)
counts<-counts[-1]
counts<-counts[which(counts>=34)]
ibol<-ibol[ibol$order_name %in% names(counts),]
ibol$order_name<-as.character(ibol$order_name)
counts<-table(ibol$order_name)
samps<-strata(ibol,"order_name",size=ceiling(0.15*(counts)))

n<-nrow(ibol)
PlayDataset<-ibol[sample.int(n,5000),]
rm(ibol)
gc()


PlayDataset<-read.csv("ibol_Play.csv")

dna<-PlayDataset$nucleotides
names(dna)<-rownames(PlayDataset)
library("Biostrings")
ibolSeq1<-DNAStringSet(PlayDataset$nucleotides)

writeXStringSet(ibolSeq1, "ibol_seqs1.fasta", append=F, format="fasta")
#write.csv(PlayDataset,file="ibol_Play.csv")

Kmer1<-oligonucleotideFrequency(ibolSeq1,1)
Kmer2<-oligonucleotideFrequency(ibolSeq1,2)
Kmer3<-oligonucleotideFrequency(ibolSeq1,3)
Kmer4<-oligonucleotideFrequency(ibolSeq1,4)
Kmer5<-oligonucleotideFrequency(ibolSeq1,5)
Kmer6<-oligonucleotideFrequency(ibolSeq1,6)
Kmer7<-oligonucleotideFrequency(ibolSeq1,7)
#Kmer8<-oligonucleotideFrequency(ibolSeq1,8)
#Kmer9<-oligonucleotideFrequency(ibolSeq1,9)
Kmers<-cbind(Kmer1,Kmer2,Kmer3,Kmer4,Kmer5,Kmer6,Kmer7)


library("energy")

O<-factor(PlayDataset$order_name)
order_dums<-model.matrix(~O-1)
norder<-ncol(order_dums)
F<-factor(PlayDataset$family_name)
family_dums<-model.matrix(~F-1)
nfamily<-ncol(family_dums)
S<-factor(PlayDataset$subfamily_name)
subfamily_dums<-model.matrix(~S-1)
nsubfamily<-ncol(subfamily_dums)
G<-factor(PlayDataset$genus_name)
genus_dums<-model.matrix(~G-1)
ngenus<-ncol(genus_dums)

dums<-cbind(order_dums,family_dums,genus_dums)

correl<-function(x,mat){
  apply(mat,2,cor,y=x)
}

dcorrel<-function(x,mat){
  apply(mat,2,dcor,y=x)
}

test<-apply(Kmers,2,correl,mat=dums)
dtest<-apply(Kmers,2,dcorrel,mat=dums)

ibolKmer11<-oligonucleotideFrequency(ibolSeq[1],11)
ibolKmer10<-oligonucleotideFrequency(ibolSeq[1],10)
ibolKmer9<-oligonucleotideFrequency(ibolSeq[1],9)
ibolKmer8<-oligonucleotideFrequency(ibolSeq[1],8)
ibolKmer7<-oligonucleotideFrequency(ibolSeq[1],7)
ibolKmer6<-oligonucleotideFrequency(ibolSeq,6)
ibolKmer5<-oligonucleotideFrequency(ibolSeq,5)
ibolKmer4<-oligonucleotideFrequency(ibolSeq,4)
ibolKmer3<-oligonucleotideFrequency(ibolSeq,3)
ibolKmer2<-oligonucleotideFrequency(ibolSeq,2)
ibolKmer1<-oligonucleotideFrequency(ibolSeq,1)
ibolKmers<-cbind(ibolKmer1,ibolKmer2,ibolKmer3,ibolKmer4,ibolKmer5,ibolKmer6)

specinf<-as.data.frame(as.character(ibol[,c("species_name","genus_name","subfamily_name","family_name","order_name","country")]))


ibolFeatures<-cbind(specinf,ibolKmers)
i <- sapply(ibolFeatures, is.factor)
ibolFeatures[i] <- lapply(ibolFeatures[i], as.character)

ibolFeatures$subfamily_name[is.na(ibol$subfamily_taxID)]<-"?"
names<-ibolFeatures[,2:5]
hierarchy<-paste(names[,1],names[,3],names[,4],sep="/")
ibolFeatures2<-cbind(ibolFeatures$species_name,hierarchy,ibolFeatures[,c(-1:-6)])

setwd("D:/Users/Dinnage")
load("ibol_feat.RData")
options( java.parameters = "-Xmx32g" )
library( "RWeka" )

write.arff(ibolFeatures2,file="ibol.arff")
write.csv(ibolFeatures2,file="ibol_Feat.csv")

save(ibolFeatures2,file="ibol_feat.RData")
save(ibolKmer1,ibolKmer2,ibolKmer3,ibolKmer4,ibolKmer5,ibolKmer6,file="ibolKmers.RData")
load("ibolKmers.RData")

get_one<-function(x){
  gc()
  x[sample.int(nrow(x),1),]
}

cl<-makeCluster(4,"SOCK")
registerDoSNOW(cl)

GeneraSet<-ddply(ibol,~genus_name,get_one,.progress="text",.parallel=F)

ibol$nucleotides<-gsub("-","",ibol$nucleotides,fixed=T)

library(data.table)
library(nnet)
ibol2<-data.table(ibol)

#GeneraSetrows<-ibol2[,{sample(.I,1)},by=genus_name]
#GeneraSet<-ibol2[GeneraSetrows[[2]]]

GeneraSetrows<-ibol2[,{.I[which.is.max(nchar(nucleotides))]},by=genus_name]
GeneraSet<-ibol2[GeneraSetrows[[2]]]


#SpeciesSetrows<-ibol2[,{sample(.I,1)},by=species_name]
#SpeciesSet<-ibol2[SpeciesSetrows[[2]]]

SpeciesSetrows<-ibol2[,{.I[which.is.max(nchar(nucleotides))]},by=species_name]
SpeciesSet<-ibol2[SpeciesSetrows[[2]]]

stopCluster(cl)

SpeciesSet2<-as.data.frame(SpeciesSet)
setwd("D:/Users/Dinnage")
write.csv(SpeciesSet2,"ibol_species1.csv")
SpeciesSet2<-read.csv("ibol_species1.csv")

SpeciesSet2$nucleotides<-gsub("-","",SpeciesSet2$nucleotides,fixed=T)
SpeciesSet2$order_name<-as.character(SpeciesSet2$order_name)
counts<-table(SpeciesSet2$order_name)
counts<-counts[-1]
counts<-counts[which(counts>=34)]
SpeciesSet<-SpeciesSet2[SpeciesSet2$order_name %in% names(counts),]
SpeciesSet<-SpeciesSet[order(SpeciesSet$order_name),]
counts<-table(SpeciesSet$order_name)

samps<-strata(SpeciesSet,"order_name",ceiling(counts*0.15))
samples<-getdata(SpeciesSet,samps)

write.csv(samples,file="COI_Play_Dataset.csv")

dna<-SpeciesSet2$nucleotides
names(dna)<-SpeciesSet2$processid
library("Biostrings")
ibolSeq1<-DNAStringSet(dna)

writeXStringSet(ibolSeq1, "ibol_seqs1.fasta", append=F, format="fasta")
#write.csv(PlayDataset,file="ibol_Play.csv")

Kmer1<-oligonucleotideFrequency(ibolSeq1,1)
Kmer2<-oligonucleotideFrequency(ibolSeq1,2)
Kmer3<-oligonucleotideFrequency(ibolSeq1,3)
Kmer4<-oligonucleotideFrequency(ibolSeq1,4)
Kmer5<-oligonucleotideFrequency(ibolSeq1,5)
Kmer6<-oligonucleotideFrequency(ibolSeq1,6)
Kmer7<-oligonucleotideFrequency(ibolSeq1,7)
#Kmer8<-oligonucleotideFrequency(ibolSeq1,8)
#Kmer9<-oligonucleotideFrequency(ibolSeq1,9)
Kmers<-cbind(Kmer1,Kmer2,Kmer3,Kmer4,Kmer5,Kmer6,Kmer7)

rm(Kmer1,Kmer2,Kmer3,Kmer4,Kmer5,Kmer6,Kmer7,Kmers)
gc()

library("energy")

O<-factor(SpeciesSet2$order_name)
order_dums<-model.matrix(~O-1)
norder<-ncol(order_dums)
F<-factor(SpeciesSet2$family_name)
family_dums<-model.matrix(~F-1)
nfamily<-ncol(family_dums)
S<-factor(SpeciesSet2$subfamily_name)
subfamily_dums<-model.matrix(~S-1)
nsubfamily<-ncol(subfamily_dums)
G<-factor(SpeciesSet2$genus_name)
genus_dums<-model.matrix(~G-1)
ngenus<-ncol(genus_dums)

dums<-cbind(order_dums,family_dums,genus_dums)

correl<-function(x,mat){
  apply(mat,2,cor,y=x)
}

dcorrel<-function(x,mat){
  apply(mat,2,dcor,y=x)
}

colsums<-apply(Kmers,2,sum)
tt1<-dcorrel(Kmers[,21000],dums)

test<-apply(Kmers,2,correl,mat=dums)
dtest<-apply(Kmers,2,dcorrel,mat=dums)

test<-RRF(Kmers,O,ntree=100)

library(mRMRe)

testd<-mRMR.data(data = as.data.frame(cbind(O,Kmers)))

data<-cbind(SpeciesSet2,Kmers)


rm(Kmers,SpeciesSet2,dums,SpeciesSet,ibol,ibol2)
gc()


#Return Time Distributions

get_RTs<-function(x,n){
  x1<-x[2:length(x)]
  x2<-x[1:(length(x)-1)]
  newx<-(x1-x2-n)
  return(newx)
}
get_RTD<-function(x){
  return(c(mean(x),var(x)))
}

get_RTDs<-function(x,data){
  n<-nchar(x)
  try<-gregexpr(x,data,fixed=T)
  test<-sapply(try,get_RTs,n=n)
  RTDs<-sapply(test,get_RTD)
  return(RTDs)
}
try<-gregexpr(kmernames[10],ibolSeq1,fixed=T)
test<-sapply(try,get_RTs,n=nchar(kmernames[10]))

RTDs<-sapply(test,get_RTD)

kmernames<-c(colnames(Kmer1),colnames(Kmer2),colnames(Kmer3),colnames(Kmer4),colnames(Kmer5),colnames(Kmer6))
RTDs<-ldply(kmernames,get_RTDs,data=ibolSeq1,.progress="text")
RTDnames<-dat[,c("processid","order_name","family_name","subfamily_name","genus_name","species_name")]

save(RTDs2,file="RTDs.RData")
save(RTDnames,fil)

setwd("D:/Users/Dinnage")
load("RTDs.RData")
RTDs<-RTDs[,1:38713]
RTDs[is.na(RTDs)]<-0

RTDs2<-data.table(RTDs)
RTDs2[is.na(RTDs2)]<-0

RTDs2<-as.matrix(RTDs)
RTDs2[is.na(RTDs2)]<-0

RTD<-t(RTDs2[apply(RTDs2,1,sum)>0,])
RTD2<-t(scale(RTD))
save(RTD2,file="RTD2.rData")

setwd("D:/Users/Dinnage")
load("RTD2.RData")

RTD<-t(RTD2)
rm(RTD2)
gc()
RTDdist<-rdist(RTD)
save(RTDdist,file="RTDdist.rData")
load("RTDdist.rData")

RTDdist1<-as.dist(RTDdist)
rm(RTDdist)
gc()

bigclust<-hclust(RTDdist1,method="average")
save(bigclust,file="BigClust.rData")

treetest<-as.phylo(bigclust)
tree2<-treetest
tree2$edge.length[which.max(treetest$edge.length)]<-0
tree2$edge.length[which.max(tree2$edge.length)]<-1
tree3<-unroot(tree2)
fun<-extract.clade(treetest,60000)

SpeciesSet2<-read.csv("ibol_species1.csv")
Grouper<-SpeciesSet2[,c("order_name","family_name","subfamily_name","genus_name","species_name")]
save(Grouper,file="Grouper.rData")
rm(SpeciesSet2)
gc()

load("Grouper.rData")
hierlab<-paste(Grouper$order_name,Grouper$family_name,Grouper$genus_name,sep="/")

newdat<-cbind(RTD,hierlab)
save(newdat,file="newdat.rData")

genustest<-as.character(sample(levels(Grouper$genus_name),200))
genustestrows<-which(Grouper$genus_name %in% genustest)
famtest<-as.character(sample(levels(Grouper$family_name),25))
famtestrows<-which(Grouper$family_name %in% famtest)
speciestestrows<-sample.int(38713,2000)
testset1<-c(genustestrows,famtestrows,speciestestrows)
testset<-testset1[!duplicated(testset1)]

testdat<-newdat[testset,]
newdat2<-newdat[-testset,]

valrows<-sample.int(34567,10000)
valdat<-newdat2[valrows,]
traindat<-newdat2[-valrows,]

rm(newdat,newdat2)
gc()
save(testset,valrows,testdat,valdat,traindat,file="ibolSets.rData")

setwd("D:/Users/Dinnage")
load("ibolSets.rData")

classes<-unique(c(testdat[,ncol(testdat)],valdat[,ncol(valdat)],traindat[,ncol(traindat)]))

fline<-"@RELATION ibolDataTest"
testfeats<-c(1:(ncol(testdat)-1))
lines<-paste("@ATTRIBUTE FEAT",testfeats," numeric",sep="")
lline<-paste("@ATTRIBUTE TAXONOMY hierarchical",paste(classes,collapse=","),sep=" ")
testfile<-c(fline,"",lines,lline,"","@DATA")

write(testfile,file="ibolTest.arff",sep="\n")
write.table(testdat,file="ibolTest.arff",append=T,quote=F,row.names=F,col.names=F,sep=",")


fline<-"@RELATION ibolDataValidation"
valfeats<-c(1:(ncol(valdat)-1))
lines<-paste("@ATTRIBUTE FEAT",valfeats," numeric",sep="")
lline<-paste("@ATTRIBUTE TAXONOMY hierarchical",paste(classes,collapse=","),sep=" ")
valfile<-c(fline,"",lines,lline,"","@DATA")

write(valfile,file="ibolVal.arff",sep="\n")
write.table(valdat,file="ibolVal.arff",append=T,quote=F,row.names=F,col.names=F,sep=",")


fline<-"@RELATION ibolDataTraining"
trainfeats<-c(1:(ncol(traindat)-1))
lines<-paste("@ATTRIBUTE FEAT",trainfeats," numeric",sep="")
lline<-paste("@ATTRIBUTE TAXONOMY hierarchical",paste(classes,collapse=","),sep=" ")
trainfile<-c(fline,"",lines,lline,"","@DATA")

write(trainfile,file="ibolTrain.arff",sep="\n")
write.table(traindat,file="ibolTrain.arff",append=T,quote=F,row.names=F,col.names=F,sep=",")


#Do GRRF
library(RRF)

cl<-makeCluster(4,"SOCK")
registerDoSNOW(cl)
#parallel test
subs<-sample.int(9179,100)
testdat<-RTD[,subs]

counts<-table(Grouper$order_name)
sampcounts<-ceiling(0.62*counts)
sampcounts[counts<10]<-counts[counts<10]
rf <- foreach(ntree=rep(125, 4), .packages='RRF') %dopar%
  RRF(RTD,as.factor(Grouper$order_name), flagReg = 0) 

stopCluster(cl)

#Reduce Dataset
counts<-table(Grouper$order_name)
keep<-Grouper$order_name %in% names(counts)[!counts<34]
newRTD<-RTD[keep,]
newGrouper<-as.character(Grouper$order_name[keep])
counts<-table(newGrouper)

Family<-as.character(Grouper$family_name[keep])

ups<-names(counts)[which(counts<500)]
upcount<-counts[which(counts<500)]
for (i in 1:length(ups)){
  rat<-ceiling(500/upcount[i])
  rows<-which(newGrouper==ups[i])
  rows2<-rep(rows,rat)
  newRTD<-rbind(newRTD,newRTD[rows2,])
  newGrouper<-c(newGrouper,newGrouper[rows2])
  print(i)
}

sampcounts<-rep(500,length(counts))
#sampcounts[counts<50]<-ceiling(0.67*counts[counts<50])
tot<-sum(sampcounts)
weights<-(tot-sampcounts)/sampcounts
weights<-weights/sum(weights)
weights2<-sampcounts/(tot-sampcounts)

counts<-table(newGrouper)
tot<-sum(counts)
weights<-(tot-counts)/counts
weights<-weights/sum(weights)
cweights<-counts/tot

#test new RRF
varUsedAll<-rep(0,ncol(newRTD))
varUsedAll[1:4]<-1
rftest <- RRF(newRTD,as.factor(newGrouper), strata=newGrouper,sampsize=sampcounts, flagReg = 1, do.trace=T, ntree=30,varUsedAll=varUsedAll)

rf1 <- RRF(newRTD,as.factor(newGrouper), strata=newGrouper,sampsize=rep(min(counts),length(counts)), flagReg = 0, do.trace=T, ntree=500)

rf2 <- RRF(newRTD,as.factor(newGrouper), strata=newGrouper,sampsize=sampcounts, flagReg = 0, do.trace=T, ntree=500)
save(rf2,file="RFtest2.rData")
#rf <- RRF(testdat,as.factor(Grouper$order_name), flagReg = 0, do.trace=T,ntree=50)
load("RFtest.rData")

#splitter
testGroup<-newGrouper
testGroup[newGrouper=="Lepidoptera"]<-"Leps"
testGroup[newGrouper!="Lepidoptera"]<-"Non-Leps"
trf1 <- RRF(newRTD,as.factor(testGroup), flagReg = 0, do.trace=T, ntree=50)


impRF2 <- rf2$importance 
impRF2 <- impRF2[,"MeanDecreaseGini"] # get the importance score 
imp2 <- impRF2/(max(impRF2)) #normalize the importance scores into [0,1]
gamma <- 0.15   #A larger gamma often leads to fewer features. But, the quality of the features selected is quite stable for GRRF, i.e., different gammas can have similar accuracy performance (the accuracy of an ordinary RF using the feature subsets). See the paper for details. 
coefReg <- (1-gamma) + gamma*imp2   # each variable has a coefficient, which depends on the importance score from the ordinary RF and the parameter: gamma
grrf2 <- RRF(newRTD,as.factor(newGrouper), strata=newGrouper, sampsize=sampcounts, flagReg=1, coefReg=coefReg, do.trace=T,ntree=500)
imp2 <- grrf2$importance
imp2 <- imp2[,"MeanDecreaseGini"]
subsetGRRF2 <- which(imp2>0) # produce the indices of the features selected by GRRF
save(grrf2,file="GRFtest2.rData")
load("GRFtest2.rData")

imporder<-order(imp2,decreasing=T)
imps<-imporder[imp2[imporder]>0]

kmermean<-paste(kmernames,"_mean",sep="")
kmervar<-paste(kmernames,"_var",sep="")
kmers<-c(rbind(kmermean,kmervar))

counts<-table(newGrouper)
tot<-sum(counts)
weights<-(tot-counts)/counts
funfun<-RRF(newRTD[,subsetGRRF2],as.factor(newGrouper), flagReg = 0, do.trace=T, ntree=500)

centroids <- classDist(newRTD, newGrouper, pca=T)

#Family level models
setwd("D:/Users/Dinnage")
load("RTD2.RData")

RTD<-t(RTD2)
rm(RTD2)
gc()

library(RRF)

load("Grouper.rData")

counts<-table(Grouper$order_name)
keep<-Grouper$order_name %in% names(counts)[!counts<34]
newRTD<-RTD[keep,]
newGrouper<-as.character(Grouper$order_name[keep])
counts<-table(newGrouper)

Family<-as.character(Grouper$family_name[keep])

load("GRFtest2.rData")
imp2 <- grrf2$importance
imp2 <- imp2[,"MeanDecreaseGini"]
subsetGRRF2 <- which(imp2>0) # produce the indices of the features selected by GRRF

#omit Leps
famOrder<-newGrouper[-which(newGrouper=="Lepidoptera")]
famFamily<-Family[-which(newGrouper=="Lepidoptera")]
famRTD<-newRTD[-which(newGrouper=="Lepidoptera"),]

#Blattaria
newsubset<-subsetGRRF2
ocounts<-table(famOrder)
oind<-which(famOrder==names(ocounts)[1])
famgroup<-famFamily[oind]
oind<-oind[-which(famgroup==" ")]
famgroup<-famFamily[oind]
newcounts<-table(famgroup)
oind<-which(famgroup %in% names(newcounts)[newcounts>1])
famgroup<-famgroup[oind]
newdat<-famRTD[oind,]
newcounts<-table(famgroup)

upper<-max(newcounts)
ups<-names(newcounts)[which(newcounts!=upper)]
upcount<-newcounts[ups]
for (i in 1:length(ups)){
  rat<-ceiling(upper/upcount[i])
  rows<-which(famgroup==ups[i])
  rows2<-rep(rows,rat)
  newdat<-rbind(newdat,newdat[rows2,])
  famgroup<-c(famgroup,famgroup[rows2])
  print(i)
}
sampler<-ceiling(mean(newcounts))
sampcounts<-rep(sampler,length(newcounts))



Famrf <- RRF(newdat,as.factor(famgroup), strata=famgroup,sampsize=sampcounts, flagReg = 0, do.trace=T, ntree=1000)
varUsedAll<-rep(0,ncol(newdat))
varUsedAll[newsubset]<-1
impRF <- Famrf$importance 
impRF <- impRF[,"MeanDecreaseGini"] # get the importance score 
imp <- impRF/(max(impRF)) #normalize the importance scores into [0,1]
gamma <- 0.15   #A larger gamma often leads to fewer features. But, the quality of the features selected is quite stable for GRRF, i.e., different gammas can have similar accuracy performance (the accuracy of an ordinary RF using the feature subsets). See the paper for details. 
coefReg <- (1-gamma) + gamma*imp   # each variable has a coefficient, which depends on the importance score from the ordinary RF and the parameter: gamma
Famgrrf <- RRF(newdat,as.factor(famgroup), strata=famgroup,sampsize=sampcounts, coefReg=coefReg,flagReg = 1, do.trace=T, ntree=1000,varUsedAll=varUsedAll)
famimp <- Famgrrf$importance
famimp <- famimp[,"MeanDecreaseGini"]
famsubset <- which(famimp>0)
newsubset<-c(newsubset,famsubset[!famsubset %in% newsubset])

finalrf<-RRF(newdat[,famsubset],as.factor(famgroup), strata=famgroup,sampsize=sampcounts, flagReg = 0, do.trace=T, ntree=2000)

library(data.table)
iboldata<-data.table(cbind(Grouper,RTD))
setkey(iboldata,order_name,family_name,subfamily_name,genus_name,species_name)
ordercounts<-iboldata[,.N,by=order_name]
goodorders<-ordercounts[N>2,order_name]
iboldata<-iboldata[as.character(goodorders),]
Grouper<-as.data.frame(iboldata[,1:5,with=F])

ordermeans<-iboldata[,{lapply(.SD,mean)},by=order_name]
ordervars<-iboldata[,{lapply(.SD,var)},by=order_name]

order_mean<-as.data.frame(ordermeans[,6:ncol(ordermeans),with=F])
rownames(order_mean)<-ordermeans[,order_name]
order_vars<-as.data.frame(ordervars[,6:ncol(ordervars),with=F])
order_means<-scale(order_mean,scale=F)
order_vars_mean<-apply(order_vars,2,mean)
s_order_means2<-order_means/order_vars
s_order_means1<-t(t(order_means)/order_vars_mean)

orderkeep<-100
familykeep<-100
genuskeep<-100

orderpca1<-prcomp(s_order_means1,F,.scale=T)
importance1<-t(t(abs(orderpca1$rotation))*orderpca1$sdev)
importance<-apply(importance1,1,mean)

#orderpca2<-prcomp(s_order_means2,F)
#importance12<-t(t(abs(orderpca2$rotation))*orderpca2$sdev)
#importance2<-apply(importance12,1,mean)

famcounts<-iboldata[,.N,by=family_name]

iboldatafam<-iboldata[,{order_mean[as.character(order_name),]}]
iboldatafam<-as.matrix(iboldata[,6:ncol(iboldata),with=F])-as.matrix(iboldatafam)
iboldatafam<-data.table(cbind(Grouper,iboldatafam))

fammeans<-iboldatafam[,{lapply(.SD,mean)},by=family_name]
famvars<-iboldatafam[,{lapply(.SD,var)},by=family_name]

fam_mean<-as.data.frame(fammeans[,6:ncol(fammeans),with=F])
rownames(fam_mean)<-fammeans[,family_name]
fam_vars<-as.data.frame(famvars[,6:ncol(famvars),with=F])
fam_means<-scale(fam_mean,scale=F)
fam_vars_mean<-apply(fam_vars,2,mean,na.rm=T)
s_fam_means2<-fam_means/fam_vars
fam_vars_mean[fam_vars_mean==0]<-min(fam_vars_mean[fam_vars_mean>0])
s_fam_means1<-t(t(fam_means)/sqrt(fam_vars_mean))

fampca1<-prcomp(s_fam_means1,F,scale.=T)
famimportance1<-t(t(abs(fampca1$rotation))*fampca1$sdev)
famimportance<-apply(famimportance1,1,mean)

iboldatagen<-iboldatafam[,{fam_mean[as.character(family_name),]}]
iboldatagen<-as.matrix(iboldatafam[,6:ncol(iboldatafam),with=F])-as.matrix(iboldatagen)
iboldatagen<-data.table(cbind(Grouper,iboldatagen))

genmeans<-iboldatagen[,{lapply(.SD,mean)},by=genus_name]
genvars<-iboldatagen[,{lapply(.SD,var)},by=genus_name]

gen_mean<-as.data.frame(genmeans[,6:ncol(genmeans),with=F])
rownames(gen_mean)<-genmeans[,genus_name]
gen_vars<-as.data.frame(genvars[,6:ncol(genvars),with=F])
gen_means<-scale(gen_mean,scale=F)
gen_vars_mean<-apply(gen_vars,2,mean,na.rm=T)
s_gen_means2<-gen_means/gen_vars
gen_vars_mean[gen_vars_mean==0]<-min(gen_vars_mean[gen_vars_mean>0])
s_gen_means1<-t(t(gen_means)/sqrt(gen_vars_mean))
gen_throwaway<-apply(s_gen_means1,2,var)==0
s_gen_means1<-s_gen_means1[,!gen_throwaway]

genpca1<-prcomp(s_gen_means1,F,scale.=T)
genimportance1<-t(t(abs(genpca1$rotation))*genpca1$sdev)
genimportance<-apply(genimportance1,1,mean)

save(orderpca1,fampca1,genpca1,iboldata,iboldatafam,iboldatagen,file="Useful.rData")
