setwd("/home/din02g/Kmer Hierarchical Classification")

library("plyr")
library("doSNOW")
library("fields")
library("RRF")

RTD<-read.csv(file="Play_RTDs.csv")
Grouper<-read.csv(file="Play_Groups.csv",stringsAsFactors=F)

#make all datasets to be run through random forest
#Orders

cl <- makeCluster(4, "SOCK")
registerDoSNOW(cl)

OGroup<-Grouper$order_name

Ocounts<-table(OGroup) #count number of observations in each order
#Ocounts<-Ocounts[Ocounts>=2] #Throw out orders with only one observation
#RTD<-RTD[OGroup %in% names(Ocounts),] #Throw out observations for excluded orders
#OGroup<-OGroup[OGroup %in% names(Ocounts)]

#Ocounts<-table(OGroup) #recount number of observations in each order

RTDscale<-scale(RTD) #scale RTDs
RTDscale[is.na(RTDscale)]<-0
RTDscale[is.nan(RTDscale)]<-0

Omeans<-ddply(cbind(as.data.frame(OGroup),as.data.frame(RTDscale)),.(OGroup),colwise(mean),.parallel=T,na.rm=T) #generate means for all RTDs by order
rownames(Omeans)<-Omeans$OGroup
Omeans<-Omeans[,-1]

stopCluster(cl)

Splitter<-function(Ocounts,Omeans,gname=0,fname=0){
  splits<-list()
  alt<-NULL
  groups<-Ocounts
  gnames<-c((fname+1):(fname+1+length(groups)*2))
  taken<-rep(FALSE,length(gnames))
  todo<-NULL
  for (i in 1:(length(groups)-1)){
    fsplit<-which.max(groups) #split off largest order
    group1<-groups[fsplit]
    group2<-groups[-fsplit]
    means1<-Omeans[names(group1),]
    means2<-Omeans[names(group2),]
    dists<-rdist(means1,means2) #generate distance between largest order and others based on RTD means
    #Add most similar orders until group1 is bigger than group2
    while (sum(group1)<sum(group2)){
      nsplit<-which.min(dists)
      group1<-c(group1,group2[nsplit])
      group2<-group2[-nsplit]
      dists<-dists[-nsplit]
    }
    #take off last added order if that grouping is more even
    if (abs(sum(group1)-sum(group2))>abs(sum(group1[-length(group1)])-sum(c(group2,group1[length(group1)])))) {
      tmp<-group1[length(group1)]
      group1<-group1[-length(group1)]
      group2<-c(group2,tmp)
    }
    splits[[i]]<-list(group1=group1,group2=group2) #save groups
    names(splits)[i]<-gname
    #taken[gname]<-TRUE
    gname1<-gnames[!taken][1]
    gname2<-gnames[!taken][2]
    taken[!taken][1:2]<-c(TRUE,TRUE)
    names(splits[[i]])<-c(gname1,gname2)
    todo<-todo[todo!=gname]
    todo<-c(todo,gname1,gname2)
    #Move to left or right group based on whether there is multiple orders in them, then repeat splitting, when get to groups of size one move to other branch
    if (length(group1)==1){
      todo<-todo[todo!=gname1]      
    }
    if (length(group2)==1){
      todo<-todo[todo!=gname2] 
    }
    gname<-todo[1]
    groups<-unlist(splits,F)[sapply(strsplit(names(unlist(splits,F)),".",T),function(x) x[2])==gname][[1]]
    #print(i)
  }
  return(splits)
}

Osplits<-Splitter(Ocounts,Omeans)
lname<-as.numeric(names(Osplits[[length(Osplits)]])[2])




#Make Order datasets

DataSplit<-function(x,gnames,dat,grouper){
  g1<-names(x[[1]])
  g2<-names(x[[2]])
  tgroups<-c(g1,g2)
  subsetter<-grouper %in% tgroups
  subdat<-dat[subsetter,]
  newgroup<-grouper[subsetter]
  newgroup[newgroup %in% g1]<-gnames[1]
  newgroup[newgroup %in% g2]<-gnames[2]
  newdat<-cbind(newgroup,subdat)
  colnames(newdat)[1]<-"OGroup"
  return(newdat)
}

DataSets<-mapply(DataSplit,Osplits,lapply(Osplits,names),MoreArgs=list(dat=RTD,grouper=OGroup),SIMPLIFY=F)

pairs<-NULL
npairs<-matrix(nrow=2,ncol=2)
code<-matrix(nrow=length(Ocounts),ncol=2)
pos<-1
for (i in 1:(length(Osplits))){
  npairs[1,]<-c(names(Osplits)[[i]],names(Osplits[[i]])[1])
  npairs[2,]<-c(names(Osplits)[[i]],names(Osplits[[i]])[2])
  for (j in 1:2){
    if (length(Osplits[[i]][[j]])==1){
      code[pos,]<-c(names(Osplits[[i]])[j],names(Osplits[[i]][[j]]))
      pos<-pos+1
    }
  }
  pairs<-rbind(pairs,npairs)
}

probFind<-list()
for (i in 1:nrow(code)){
  con=TRUE
  spec<-code[i,1]
  rs<-NULL
  while (con==TRUE){
    tt<-which(pairs[,2]==spec)
    if(any(tt)){
      rr<-pairs[tt,]
      rs<-c(rs,tt)
      spec<-rr[1]
    }else{con<-FALSE}
  }
  probFind[[i]]<-pairs[rs,]
}
names(probFind)<-code[,2]

#Family
lname<-as.numeric(names(Osplits[[length(Osplits)]])[2])
FDataSets<-list()
Fsplits<-list()
Fpairs<-list()
Fcode<-list()
FprobFind<-list()

GDataSets<-list()
Gsplits<-list()
Gpairs<-list()
Gcode<-list()
GprobFind<-list()

for (i in 1:length(Ocounts)){
  Ord<-names(Ocounts)[i]
  Fdat<-RTD[Grouper$order_name==Ord,]
  FGrouper<-Grouper[Grouper$order_name==Ord,]
  FGroup<-FGrouper$family_name
  Fcounts<-table(FGroup) #count number of observations in each family
  Fcounts<-Fcounts[Fcounts>=2] #Throw out families with only one observation
  Fdat<-Fdat[FGroup %in% names(Fcounts),] #Throw out observations for excluded families
  FGroup<-FGroup[FGroup %in% names(Fcounts)]
  Fcounts<-table(FGroup) #recount number of observations in each family
  if (length(Fcounts)>1){
    Fdatscale<-scale(Fdat) #scale RTDs
    Fdatscale[is.na(Fdatscale)]<-0
    Fdatscale[is.nan(Fdatscale)]<-0
    cl <- makeCluster(4, "SOCK")
    registerDoSNOW(cl)
    Fmeans<-ddply(cbind(as.data.frame(FGroup),as.data.frame(Fdatscale)),.(FGroup),colwise(mean),.parallel=T,na.rm=T) #generate means for all RTDs by family
    stopCluster(cl)
    rownames(Fmeans)<-Fmeans$FGroup
    Fmeans<-Fmeans[,-1]
    Oname<-as.numeric(code[code[,2]==Ord,1])
    Splits<-Splitter(Fcounts,Fmeans,Oname,lname)
    lname<-as.numeric(names(Splits[[length(Splits)]])[2])
    
    FDataSets[[Ord]]<-mapply(DataSplit,Splits,lapply(Splits,names),MoreArgs=list(dat=Fdat,grouper=FGroup),SIMPLIFY=F)
    
    Pairs<-NULL
    npairs<-matrix(nrow=2,ncol=2)
    Code<-matrix(nrow=length(Fcounts),ncol=2)
    pos<-1
    for (i in 1:(length(Splits))){
      npairs[1,]<-c(names(Splits)[[i]],names(Splits[[i]])[1])
      npairs[2,]<-c(names(Splits)[[i]],names(Splits[[i]])[2])
      for (j in 1:2){
        if (length(Splits[[i]][[j]])==1){
          Code[pos,]<-c(names(Splits[[i]])[j],names(Splits[[i]][[j]]))
          pos<-pos+1
        }
      }
      Pairs<-rbind(Pairs,npairs)
    }
    
    ProbFind<-list()
    for (i in 1:nrow(Code)){
      con=TRUE
      spec<-Code[i,1]
      rs<-NULL
      while (con==TRUE){
        tt<-which(Pairs[,2]==spec)
        if(any(tt)){
          rr<-Pairs[tt,]
          rs<-c(rs,tt)
          spec<-rr[1]
        }else{con<-FALSE}
      }
      ProbFind[[i]]<-Pairs[rs,]
    }
    names(ProbFind)<-Code[,2]
    Fsplits[[Ord]]<-Splits
    Fpairs[[Ord]]<-Pairs
    Fcode[[Ord]]<-Code
    FprobFind[[Ord]]<-ProbFind
    
    GDataSets1<-list()
    Gsplits1<-list()
    Gpairs1<-list()
    Gcode1<-list()
    GprobFind1<-list()
    
    for (i in 1:length(Fcounts)){
      Fam<-names(Fcounts)[i]
      Gdat<-RTD[Grouper$family_name==Fam,]
      GGrouper<-Grouper[Grouper$family_name==Fam,]
      GGroup<-GGrouper$genus_name
      Gcounts<-table(GGroup) #count number of observations in each genus
      Gcounts<-Gcounts[Gcounts>=2] #Throw out genera with only one observation
      Gdat<-Gdat[GGroup %in% names(Gcounts),] #Throw out observations for excluded genera
      GGroup<-GGroup[GGroup %in% names(Gcounts)]
      Gcounts<-table(GGroup) #recount number of observations in each genera
      if (length(Gcounts)>1){
        Gdatscale<-scale(Gdat) #scale RTDs
        Gdatscale[is.na(Gdatscale)]<-0
        Gdatscale[is.nan(Gdatscale)]<-0
        cl <- makeCluster(4, "SOCK")
        registerDoSNOW(cl)
        Gmeans<-ddply(cbind(as.data.frame(GGroup),as.data.frame(Gdatscale)),.(GGroup),colwise(mean),.parallel=T,na.rm=T) #generate means for all RTDs by genera
        stopCluster(cl)
        rownames(Gmeans)<-Gmeans$GGroup
        Gmeans<-Gmeans[,-1]
        Fname<-as.numeric(Code[Code[,2]==Fam,1])
        Splits1<-Splitter(Gcounts,Gmeans,Fname,lname)
        lname<-as.numeric(names(Splits1[[length(Splits1)]])[2])
        
        GDataSets1[[Fam]]<-mapply(DataSplit,Splits1,lapply(Splits1,names),MoreArgs=list(dat=Gdat,grouper=GGroup),SIMPLIFY=F)
        
        Pairs1<-NULL
        npairs<-matrix(nrow=2,ncol=2)
        Code1<-matrix(nrow=length(Gcounts),ncol=2)
        pos<-1
        for (i in 1:(length(Splits1))){
          npairs[1,]<-c(names(Splits1)[[i]],names(Splits1[[i]])[1])
          npairs[2,]<-c(names(Splits1)[[i]],names(Splits1[[i]])[2])
          for (j in 1:2){
            if (length(Splits1[[i]][[j]])==1){
              Code1[pos,]<-c(names(Splits1[[i]])[j],names(Splits1[[i]][[j]]))
              pos<-pos+1
            }
          }
          Pairs1<-rbind(Pairs1,npairs)
        }
        
        ProbFind1<-list()
        for (i in 1:nrow(Code1)){
          con=TRUE
          spec<-Code1[i,1]
          rs<-NULL
          while (con==TRUE){
            tt<-which(Pairs1[,2]==spec)
            if(any(tt)){
              rr<-Pairs1[tt,]
              rs<-c(rs,tt)
              spec<-rr[1]
            }else{con<-FALSE}
          }
          ProbFind1[[i]]<-Pairs1[rs,]
        }
        names(ProbFind1)<-Code1[,2]
        Gsplits1[[Fam]]<-Splits1
        Gpairs1[[Fam]]<-Pairs1
        Gcode1[[Fam]]<-Code1
        GprobFind1[[Fam]]<-ProbFind1
      }
      print(Fam)
    }
    Gsplits[[Ord]]<-Gsplits1
    Gpairs[[Ord]]<-Gpairs1
    Gcode[[Ord]]<-Gcode1
    GprobFind[[Ord]]<-GprobFind1
  }
  print(Ord)
}

FData<-do.call(c,FDataSets)
names(FData)<-sapply(strsplit(names(FData),".",T),function(x) x[2])
nDataSets<-c(DataSets,FData)



Do_RF<-function(x,ntree=500,do.trace=T,features=NULL){
  grouper<-as.factor(x[,1])
  dat<-x[,-1]
  if (!is.null(features)){
    dat<-dat[,features]
  }
  samps<-rep(min(table(grouper)),2)
  RF<-RRF(dat,grouper,ntree=ntree,strata=grouper,sampsize=samps,flagReg=0,do.trace=do.trace)
  return(RF)
}

#rfs<-llply(DataSets,Do_RF,ntree=500)
nrfs<-llply(nDataSets,Do_RF,ntree=500)

#grrfs<-list()
#feaIni<-NULL
#for (i in 1:length(rfs)){
#  grouper<-as.factor(DataSets[[i]][,1])
#  dat<-DataSets[[i]][,-1]
#  samps<-rep(min(table(grouper)),2)
#  impRF=rfs[[i]]$importance
#  impRF=impRF[,"MeanDecreaseGini"]
#  imp=impRF/(max(impRF))#normalize the importance score
#  gamma = 0.25
#  coefReg=(1-gamma)+gamma*imp #weighted average
#  grrf <- RRF(dat,grouper,coefReg=coefReg, flagReg=1, feaIni=feaIni, do.trace=T)
#  grrfs[[i]]<-grrf
#  feaIni<-grrf$feaSet
#}

#feat<-grrfs[[length(grrfs)]]$feaSet
#grfs<-llply(DataSets,Do_RF,ntree=500,features=feat)

ngrrfs<-list()
feaIni<-NULL
for (i in 1:length(nrfs)){
  grouper<-as.factor(nDataSets[[i]][,1])
  dat<-nDataSets[[i]][,-1]
  samps<-rep(min(table(grouper)),2)
  impRF=nrfs[[i]]$importance
  impRF=impRF[,"MeanDecreaseGini"]
  imp=impRF/(max(impRF))#normalize the importance score
  gamma = 0.25
  coefReg=(1-gamma)+gamma*imp #weighted average
  ngrrf <- RRF(dat,grouper,coefReg=coefReg, flagReg=1, feaIni=feaIni, do.trace=T)
  ngrrfs[[i]]<-grrf
  feaIni<-ngrrf$feaSet
}

nfeat<-ngrrfs[[length(ngrrfs)]]$feaSet
ngrfs<-llply(nDataSets,Do_RF,ntree=500,features=nfeat)

#predictions<-llply(grfs,predict,newdata=RTD[,feat],type="prob")
npredictions<-llply(ngrfs,predict,newdata=RTD[,nfeat],type="prob")

cmult<-function(x){
  mult<-x[,1]
  if (ncol(x)>1){
    for (i in 2:ncol(x)){
      mult<-mult*x[,i]
    }
  }
  return(mult)
}
predictTaxa<-function(x,pred,probFind){
    probF<-matrix(probFind[x][[1]],ncol=2)
    preds<-apply(probF,1,function(y){pred[[y[1]]][,y[2]]})
    print(x)
    return(cmult(preds))
}

Fprob<-do.call(c,FprobFind)
names(Fprob)<-sapply(strsplit(names(Fprob),".",T),function(x) x[2])
nprobFind<-c(probFind,Fprob)

#predTaxa<-ldply(names(probFind),predictTaxa,pred=predictions,probFind=probFind)
#predTaxa<-t(predTaxa)
#colnames(predTaxa)<-names(probFind)

npredTaxa<-ldply(names(nprobFind),predictTaxa,pred=npredictions,probFind=nprobFind)
npredTaxa<-t(npredTaxa)
colnames(npredTaxa)<-names(nprobFind)

FFcounts<-table(Grouper$family_name)
FFcounts<-FFcounts[FFcounts>=2]
FGrouper<-Grouper[Grouper$family_name %in% names(FFcounts),]
TT<-table(FGrouper[,c("order_name","family_name")])
Of<-apply(TT,2,function(x) rownames(TT)[which(x>0)])
cutmult<-do.call(rbind,lapply(Opredicts,function(x) Of==x))

Opred<-npredTaxa[,1:length(Ocounts)]
Ocuts<-Opred
Ocuts<-
OwhichTaxa<-apply(Ocuts,1,which.max)
Opredicts<-colnames(Opred)[OwhichTaxa]
Otrue_pred<-cbind(OGroup,Opredicts)
colnames(Otrue_pred)<-c("True_Taxa","Pred_Taxa")
Oconmat<-ddply(as.data.frame(Otrue_pred),"True_Taxa",function (x) table(x[,2]))

Fpred<-npredTaxa[,(length(Ocounts)+1):ncol(npredTaxa)]
Fcuts<-Fpred*cutmult[,colnames(Fpred)]
FwhichTaxa<-apply(Fcuts,1,which.max)
Fpredicts<-colnames(Fpred)[FwhichTaxa]
FFGroup<-Grouper$family_name
Ftrue_pred<-cbind(FFGroup,Fpredicts)
colnames(Ftrue_pred)<-c("True_Taxa","Pred_Taxa")
Fconmat<-ddply(as.data.frame(Ftrue_pred),"True_Taxa",function (x) table(x[,2]))

Fsumms<-list()
Fcon<-Fconmat[,-1]
rownames(Fcon)<-Fconmat$True_Taxa
Fcon<-Fcon[colnames(Fcon),]
for (i in 1:nrow(Fcon)){
  Fsumms[[i]]<-Fcon[i,i]/sum(Fcon[i,])
}
Fsumms<-unlist(Fsumms)

Osumms<-list()
Ocon<-conmat[,-1]
rownames(Ocon)<-conmat$True_Taxa
Ocon<-Ocon[colnames(Ocon),]
for (i in 1:nrow(Ocon)){
  Osumms[[i]]<-Ocon[i,i]/sum(Ocon[i,])
}
Osumms<-unlist(Osumms)

#Comparison
cut<-Ocounts/sum(Ocounts)
nrf<-RRF(RTD,as.factor(OGroup),cutoff=cut,flagRed=0,do.trace=T)

impRF=nrf$importance
impRF=impRF[,"MeanDecreaseGini"]
imp=impRF/(max(impRF))#normalize the importance score
gamma = 0.25
coefReg=(1-gamma)+gamma*imp #weighted average
ngrrf <- RRF(RTD,as.factor(OGroup),cutoff=cut,coefReg=coefReg, flagReg=1, do.trace=T)
feat<-ngrrf$feaSet
ngrf<-RRF(RTD[,feat],as.factor(OGroup),cutoff=cut,flagRed=0,do.trace=T)

#Comparison No Cutoff

nrf2<-RRF(RTD,as.factor(OGroup),flagRed=0,do.trace=T)

impRF=nrf2$importance
impRF=impRF[,"MeanDecreaseGini"]
imp=impRF/(max(impRF))#normalize the importance score
gamma = 0.25
coefReg=(1-gamma)+gamma*imp #weighted average
ngrrf2 <- RRF(RTD,as.factor(OGroup),coefReg=coefReg, flagReg=1, do.trace=T)
feat2<-ngrrf$feaSet
ngrf2<-RRF(RTD[,feat2],as.factor(OGroup),flagRed=0,do.trace=T)
