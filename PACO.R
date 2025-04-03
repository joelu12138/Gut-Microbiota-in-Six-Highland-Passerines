#######The following codes reference published literature
#####Youngblut, Nicholas D., et al. "Host diet and evolutionary history explain different aspects of gut microbiome diversity among vertebrate clades." Nature communications 10.1 (2019): 1-15.
#work_dir="/media/project/2.R/LIPA"
#setwd(work_dir)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ape)
library(paco)
library(picante)
library(future)
library(future.batchtools)
library(doParallel)
rescale_dist_mtx = function(m){
  m = m %>% as.matrix
  labs = m %>% colnames
  n_row = m %>% nrow
  n_col = m %>% ncol
  x = m %>% as.vector 
  x = scales::rescale(x) 
  m = matrix(x, nrow=n_row, ncol=n_col)
  colnames(m) = labs
  rownames(m) = labs
  m = m %>% as.dist
  return(m)
}
## load files
host_tree = read.tree("Phylogenetic relationships of the host.nwk")
host_tree =multi2di(host_tree)   #convert to rooted tree
#loading the relative abundance of microbial genera
spe<-read.csv("OTU_table.csv",header = T,row.names = 1)
head(spe)
spe<-t(spe)
#loading phylogenetic relationships of microbial genera
tree<-read.tree("Phylogenetic relationships of microorganisms.tree")
taxa_names(tree)
tree_tips <- tree$tip.label
species_in_df <- colnames(spe)
tips_to_drop <- setdiff(tree_tips, species_in_df)
pruned_tree <- drop.tip(tree, tips_to_drop)
tree_tips <- pruned_tree$tip.label
df_pruned <- spe[, tree_tips]
dim(df_pruned)
View(df_pruned)
write.tree(pruned_tree,"pruned_tree.tre")
write.csv(df_pruned,"MIC_tree.csv")

####
##16S dist
otu<-read.csv("MIC_tree.csv",header = T,row.names = 1)
head(otu)
otu<-t(otu)
#
phy.tree<-pruned_tree
micro_D = phy.tree %>% cophenetic %>% rescale_dist_mtx %>% as.matrix
micro_D %>% dim
##host dist
host_D = host_tree %>% cophenetic %>% rescale_dist_mtx %>% as.matrix
host_D %>% dim
###Preparing data
otu = otu %>% t %>% apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix
D = prepare_paco_data(H=host_D, P=micro_D, HP=df_pruned)
D %>% names
D = add_pcoord(D, correction='cailliez')
D %>% names
options(future.globals.maxSize = 5 * 1024^3)  # Set the limit to 5 GiB
#########
###PACO: diffused model
PACo_file = file.path('physeq_IndD_PACo-Con.RDS')
D %<-% { PACo(D, nperm=1000, seed=3874, method='quasiswap', symmetric=TRUE) } %packages% c("paco")
#
sink("res+PACO.txt")
D
sink()