require(rPython)
require(plyr)
require(dplyr)
require(magrittr)

python.load("/home/din02g/Al-Fr_Diversity/Code/SourceCode.py")
python.load("/home/din02g/Al-Fr_Diversity/Code/suffix_tree.py")
python.exec("import sys")
python.exec("sys.path.append(\"/home/din02g/Al-Fr_Diversity/Code/\")")

getRTDs<-function(path,k=5){
  rtds<-python.call("computeRTDfromfastafile",path,k) 
  out<- rtds[[2]] %>% laply(identity)
  rtdnames<-python.call("returnNames",k)
  newnames<-laply(rtdnames,function(x) paste(x,c("mean","var"),sep="_"))
  rownames(out)<-rtds[[1]]
  colnames(out)<-newnames
  return(out)
}


getAll_RTDs<-function(path, k=5) {
  out<-llply(seq_len(k),function(x) getRTDs(path,x)) %>% do.call(cbind,.)
  return(out)
}

