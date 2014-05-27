
#dat<-read.csv("D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
#dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
#write.csv(dat,file="D:/Users/Dinnage/Projects/WBHC-Project/data/COI_Play_Dataset_clean.csv")


dat<-read.csv("/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset.csv",stringsAsFactors=F)
dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
write.csv(dat,file="/home/din02g/Google Drive/WBHC-Project/data/COI_Play_Dataset_clean.csv")

dat<-read.csv("/home/din02g/Google Drive/WBHC-Project/data/Train/COI_Train_Dataset.csv",stringsAsFactors=F)
dat<-dat[,c("processid","nucleotides","order_name","family_name","subfamily_name","genus_name","species_name")]
write.csv(dat,file="/home/din02g/Google Drive/WBHC-Project/data/Train/COI_Train_Dataset_clean.csv")