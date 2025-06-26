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

#Leave-One-Out Kruskal-Wallis analysis
results_loo <- sapply(1:nrow(data), function(i) {
  data_loo <- data[-i, ]
  test <- kruskal_test(data_loo, Metazoa ~ Group)
  return(test$p)
})

plot(results_loo, type = "b", main = "Leave-One-Out Kruskal p-values",
     xlab = "Removed Sample Index", ylab = "p-value")
