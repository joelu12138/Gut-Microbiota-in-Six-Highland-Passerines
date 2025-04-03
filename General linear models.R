#loading package
library(tidyverse)
library(colorspace)
library(ggplot2)
#laoding OTU table
otu <-read.csv("OTU_table.csv",header = T,row.names = 1)#OTU table of diet or gut microbiota
head(otu)
otu <- t(otu)
head(otu)
otu.distance <- vegdist(otu, method = 'bray')
otu.distance
write.csv(otu.distance,"otu.distance.csv")
#loading data
mat<-read.csv("otu.distance.csv",header = T,row.names = 1)
head(mat)
long_data <- as.data.frame(mat) %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%
  filter(match(Row, rownames(mat)) > match(Column, colnames(mat)))
head(long_data)
write.csv(long_data,"data2.csv")
#
data<-read.csv("data2.csv",header = T)
head(data)
data2=scale(data)
data3=as.data.frame(data2)
write.csv(data3,"data4.csv")
head(data3)
#
data<-read.csv("data4.csv",header = T,row.names = 1)
head(data)
result<-lm(MIC~Diet,data=data)
summary(result)
fit <-lm(MIC~Diet+Taxa,data=data)
summary(fit)
#
slope=7.181e-04
intercept=0.990
pdf("CorrHt2.pdf", height = 8, width = 8)
ggplot()+
  geom_point(data=data,aes(x=Fungi,y=p__Euryarchaeota),size=5,stroke=0.3,shape=21,fill="firebrick3")+
  geom_smooth(data=data,aes(x=Fungi,y=p__Euryarchaeota),method="lm",formula="y~x",color="firebrick3",se=T,size=2)+
  geom_abline(slope = slope, intercept = intercept, color = "#293890", linetype = "dashed", size = 1.5) +  # the result of considering with host taxa
  labs(x=expression("Phylogenetic distance"),y=expression("gut microflora diversity"))+
  theme_test()+#
  theme(plot.margin=unit(c(1,1,1,1),"pt"),
        legend.position="non",
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.key.size=unit(1,"cm"),
        legend.background=element_blank(),
        axis.text=element_text(color="black",size=8),
        axis.text.x=element_text(size=8),
        axis.title=element_text(size=9),
        title=element_text(size=11))+
  annotate("text",label=expression(italic(r)*"=-6.955e-04,"*italic(P)*"=0.990"),size=10)+
  annotate("text",label=expression(italic(r)*"=7.181e-04,"*italic(P)*"=0.990"),size=10)
dev.off()