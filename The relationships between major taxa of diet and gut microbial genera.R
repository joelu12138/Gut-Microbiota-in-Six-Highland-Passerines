#loading data
MIC<-read.csv("Relative abundance of gut microbial genera.csv",header = T)#Relative abundance of gut microbial genera or phyla
Diet<-read.csv("Relative abundance of diet maior taxa.csv",header = T)
data<-merge(MIC,Diet,by="samples")
write.csv(data,"data.csv")
#
data=read.csv("data,csv",header = T,row.names = 1)
data2=scale(data)
data3=as.data.frame(data2)
write.csv(data3,"data4.csv")
data3=read.csv("data4.csv",header = T,row.names = 1)
head(data3)
dim(data3)
regression2 <- function(x){
  coef(summary(lm(x~SAR+Group,data=data3)))[2,c(1,2,3,4)]
}
result <- map(data3[,7:22],regression2) %>% 
  do.call(rbind,.) %>% 
  as.data.frame() %>% 
  rownames_to_column(.)
write.csv(result,"SAR-lm.csv")
#
#Heatmap
library(psych)
library(reshape2)
library(pheatmap)
pmt=read.csv("pvalue.csv",header = T,row.names = 1)
head(pmt)
cmt=read.csv("cor.csv",header = T,row.names = 1)
head(cmt)
#
if (!is.null(pmt)) {
  ssmt <- pmt < 0.01
  pmt[ssmt] <- '**'
  smt <- pmt > 0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt] <- ''
} else {
  pmt <- F
}
pmt
#
pheatmap(cmt, scale = "none", cluster_row = F, cluster_col = F, gaps_row = c(4),               
         display_numbers = pmt, fontsize_number = 20, number_color = "white",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(200),
         cellwidth = 20, cellheight = 20,filename="heatmap.pdf")