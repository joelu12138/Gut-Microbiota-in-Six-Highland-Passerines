##laoding packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggsci)

#laoding OTU table
asv<-read.csv("OTU_table.csv",header = T,row.names = 1)#OTU table of diet or gut microbiota
##All numbers are converted to numeric type
asv <- asv[,-ncol(asv)]
for(i in 1: ncol(asv))
{asv[,i]<-as.numeric(asv[,i])}

asv_t <- as.data.frame(t(asv))
rare <- rowSums(asv_t)
step <- 1000
alpha_rare <- list()

for (i in 1:nrow(asv_t)) {
  step_num <- seq(1, rare[i], step)
  if (max(step_num) < rare[i]){ step_num <- c(step_num, as.numeric(rare[i])) }
  alpha_rare_i <- NULL
  for (step_num_n in step_num) {
    asv_n <- rrarefy(asv_t[i, ], step_num_n)
    alpha_rare_i <- c(alpha_rare_i, diversity(asv_n, index = 'shannon', base = 2))
  }
  names(alpha_rare_i) <- step_num
  alpha_rare <- c(alpha_rare, list(alpha_rare_i))
}
names(alpha_rare) <- rownames(asv_t)
plot_data <- data.frame()
for (i in names(alpha_rare)) {
  alpha_curves_i <- (alpha_rare[[i]])
  alpha_curves_i <- data.frame(rare = names(alpha_curves_i), alpha = alpha_curves_i, sample =i, stringsAsFactors = FALSE)
  plot_data <- rbind(plot_data, alpha_curves_i)
}
plot_data$rare <- as.numeric(plot_data$rare)
plot_data$alpha <- as.numeric(plot_data$alpha)
plot_data$sample <- ordered(plot_data$sample,levels = unique(metadata$Sampleid))
#
p1 <- ggplot(plot_data, aes(rare,alpha, color = sample)) +
  geom_smooth(se = FALSE, method = "lm" ,formula = y ~ log(x)) +
  theme( panel.grid = element_blank(),
         panel.background = element_rect(fill = 'transparent', color = 'black'),legend.key =element_rect(fill = 'transparent')) +
  labs(x = "depth", y = "shannon")
ggsave(p1,filename = "Dilution curve.pdf",width = 8, height = 5)
