#loading packages
library(vegan)
library(ggplot2)
#laoding OTU table
otu <-read.csv("OTU_table.csv",header = T,row.names = 1)#OTU table of gut microbiota
head(otu)
otu <- t(otu)
head(otu)
#
otu.distance <- vegdist(otu, method = 'bray')
otu.distance
write.csv(otu.distance,"otu.distance.csv")
#calculate Dissimilarity distance
pcoa <- cmdscale (otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
group<-read.csv("Group.csv",header = T)
df <- merge(pc12,group,by="samples")
#PERMANOVA statistical test was performed and adonis2 function replacement test was used
set.seed(121314)
adonis_result<-adonis2(otu ~ Group, data = group,permutations = 9999,method="bray")
adonis_result
#
dune_adonis <- paste0("adonis R2: ",round(adonis_result$R2,4), "; P-value: ", adonis_result$`Pr(>F)`)
head(df)
p=ggplot(data=df,aes(x=V1,y=V2,color=Group))+
  theme_bw()+
  geom_point(size=8)+
  #stat_ellipse()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03, vjust=0),size=0)+
  guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PCoA1 ","(",pc[1],"%)"),
       y=paste0("PCoA2 ","(",pc[2],"%)"),
       title=dune_adonis)+
  #scale_color_manual(values = color) +
  #scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16,angle=90),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        panel.grid=element_blank())
p
ggsave(p,filename = "PCoA_bray.pdf",width = 10, height = 10)