# loading packages
library(tidyverse)
library(ggpubr)
library(rstatix)
library(PMCMRplus)
data<-read.csv("data.csv",header = T,row.names = 1)#data of ¦Á-diversity or Relative abundance of gut microbial phyla
head(data)
res.kruskal <- data %>% kruskal_test(core ~ Group)
res.kruskal
pwc2 <- data %>%
  wilcox_test(core ~ Group, p.adjust.method = "bonferroni")
pwc2
write.csv(pwc2,"Chromadoridae.csv")