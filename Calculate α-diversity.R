#
setwd("C://Rworkplace")
#
library(vegan)
library(picante)
#To read the data, header=T means that the first line in the file is set to the column name. 
#Row.names = 1 indicates that the first column is set to the row name
otu <- read.csv("OTU_table.csv",header = T,row.names = 1)#OTU table of diet or gut microbiota
head(otu)
otu <- t(otu)
## Richness
richness <- rowSums(otu > 0)
#Shannon
shannon<-diversity(otu)
#Pielou 
pielou <- shannon_index / log(richness, exp(1))
##Simpson
simpson<- diversity(otu, index = 'simpson')
#Chao1 
chao1 <- estimateR(otu)[2, ]
#ACE 
ace <- estimateR(otu)[4, ]
report=cbind(richness,shannon,pielou,simpson,chao1,ace)
head(report)
write.csv(report,"report.csv")