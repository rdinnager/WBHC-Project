setwd("/data/noflush/din02g")

dat1<-read.csv("Train_RTDs_1.csv",stringsAsFactors=F)
dat2<-read.csv("Train_RTDs_2.csv",stringsAsFactors=F)
dat3<-read.csv("Train_RTDs_3.csv",stringsAsFactors=F)
dat4<-read.csv("Train_RTDs_4.csv",stringsAsFactors=F)
dat5<-read.csv("Train_RTDs_5.csv",stringsAsFactors=F)
dat6<-read.csv("Train_RTDs_6.csv",stringsAsFactors=F)

dat<-rbind(dat1,dat2,dat3,dat4,dat5,dat6)
write.csv(dat,file="Train_RTDs.csv",row.names=F)

grp1<-read.csv("Train_Groups_1.csv",stringsAsFactors=F)
grp2<-read.csv("Train_Groups_2.csv",stringsAsFactors=F)
grp3<-read.csv("Train_Groups_3.csv",stringsAsFactors=F)
grp4<-read.csv("Train_Groups_4.csv",stringsAsFactors=F)
grp5<-read.csv("Train_Groups_5.csv",stringsAsFactors=F)
grp6<-read.csv("Train_Groups_6.csv",stringsAsFactors=F)

grp<-rbind(grp1,grp2,grp3,grp4,grp5,grp6)
write.csv(grp,file="Play_Groups.csv",row.names=F)