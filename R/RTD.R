require(rPython)
require(Biostrings)
require(plyr)
require(dplyr)
require(magrittr)

python.load("/home/din02g/Al-Fr_Diversity/Code/SourceCode.py")
python.load("/home/din02g/Al-Fr_Diversity/Code/suffix_tree.py")
python.exec("import sys")
python.exec("sys.path.append(\"/home/din02g/Al-Fr_Diversity/Code/\")")

getRTDs<-function(path,k=5){
  rtds<-python.call("computeRTDfromfastafile",path,k) 
  out<- rtds[[2]] %>% ldply(identity)
  rtdnames<-python.call("returnNames",k)
  newnames<-unlist(lapply(rtdnames,function(x) paste(x,c("mean","var"),sep="_")))
  rownames(out)<-rtds[[1]]
  colnames(out)<-newnames
  return(out)
}


getAll_RTDs<-function(path, k=5) {
  all.rtds<-llply(seq_len(k),function(x) getRTDs(path,x)) %>% 
  out<-do.call(cbind,all.rtds)
  return(out)
}

