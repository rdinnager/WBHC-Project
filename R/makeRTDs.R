#seqdat <- read.csv("data/FullData/insect_COI_data_full_July14_2014.csv")

k <- 6
fastafile <- paste(getwd(), "/data/FullData/insect_COI_data_full_July14_2014.csv", sep = "")
system(paste("python D:/Users/Dinnage/Projects/WBHC-Project/inst/scripts/RTD.py", fastafile, k, 
             "> D:/Users/Dinnage/Projects/WBHC-Project/data/FullData/insect_COI_data_full_RTDs_July17_2014.csv"))