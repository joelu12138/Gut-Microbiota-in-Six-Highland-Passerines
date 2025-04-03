#laoding data
library(vegan)
#set the working directory
setwd("C://Rworkplace")
#loading OTU table
otu<-read.csv("OTU_table.csv",header = T,row.names = 1)#OTU table of diet or gut microbiota
head(otu)
#View the sum of each sample
colSums(otu)
#Data pumping
otu_Flattening = as.data.frame(t(rrarefy(t(otu),min(colSums(otu)))))
#View the sum of each sample flattened
colSums(otu_Flattening)
#Save the flattened otu table to the working directory for later diversity analysis
write.table(otu_Flattening,file="otu_Flattening_Diet.csv",sep=",",quote=FALSE)
