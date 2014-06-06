#Script to calculate Return Time Distribution from Full Length COI Data chopped into 90 base-pair segments
library(rPython)
source("/home/din02g/Google Drive/WBHC-Project/R/RTD.R")
library(Biostrings)

## path to the segment data files
path <- "/home/din02g/Google Drive/WBHC-Project/data/FullCOI/Segments"
numfiles <- 940
maxK <- 6

## loop through all data files
for (i in 1:numfiles) {
  fullpath <- paste(path, "/Hexapoda_COI_segments_", sprintf("%03d",i),".csv", sep="")
  ## read in dataset
  dat <- read.csv(fullpath, stringsAsFactors=FALSE)
  ## pull out sequences and convert them to DNAStringSet (from Biostrings package)
  seqs <- DNAStringSet(dat$Sequence)
  ## give the sequences their ID number so that I can keep track of them later
  names(seqs) <- dat$IDnum
  ## write sequence to a fasta file for use with Python code
  fasta.name <- paste(path, "/fasta/Hexapoda_COI_segments_seqs_", sprintf("%03d",i),".fa", sep="")
  writeXStringSet(seqs, fasta.name)
  ## get RTD data using function that interfaces with Python code
  RTDs <- getAll_RTDs(fasta.name, k=maxK)
  ## attach results to original data and write to a csv for later analysis
  newdat <- cbind(dat, RTDs)
  dat.name <- paste(path, "/RTDs/Hexapoda_COI_segments_RTD_", sprintf("%03d",i),".csv", sep="")
  write.csv(newdat, file=dat.name, row.names=FALSE)
  gc()
  print(paste("Done: File number", i, "out of", numfiles))
}

## testing
seqs<-seqs[501:1000]
fasta.name <- paste(path, "/fasta/insect_COI_data_clean_seqs_", sprintf("%03d",i),".fa", sep="")
writeXStringSet(seqs, fasta.name)
system.time(
  RTDs <- getAll_RTDs(fasta.name, k=maxK)
)

## on this system there is enough memory to do 500 sequences at a time
## so I split the data into 500 row chuncks before hand, hence the above loop
