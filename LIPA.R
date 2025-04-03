#######The following codes reference published literature
#####Youngblut, Nicholas D., et al. "Host diet and evolutionary history explain different aspects of gut microbiome diversity among vertebrate clades." Nature communications 10.1 (2019): 1-15.
install.packages("phylosignal")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(phyloseq)
library(phylosignal)
library(doParallel)
#####fuctions#####
# phylosignal analysis on each tree subsample
phyloSignal_each = function(i, tree_4d, methods='all', reps=99){
  # phylosignal
  physig_res = phyloSignal(tree_4d[[i]], methods=methods, reps=reps)
  
  # formatting output
  tmp1 = physig_res$stat 
  tmp1$OTU = rownames(tmp1)
  
  # calculating qvalue on results
  tmp2 = physig_res$pvalue
  tmp3 = apply(tmp2, 2, function(x) p.adjust(x, method='BH')) %>% as.data.frame
  tmp2$OTU = rownames(tmp2)
  tmp3$OTU = rownames(tmp2)
  
  tmp1 = tmp1 %>%
    gather(method, coef, -OTU)
  tmp2 = tmp2 %>%
    gather(method, pvalue, -OTU)
  tmp3 = tmp3 %>%
    gather(method, qvalue, -OTU)
  
  tmp1 %>%
    inner_join(tmp2, c('OTU', 'method')) %>%
    inner_join(tmp3, c('OTU', 'method')) %>%
    mutate(subsample_rep = i)
}
#' randomly selecting one per group
phylo4d_subsample = function(L, df, otu, tree){
  # get subsample (one sample per species)
  df = df %>% 
    group_by(scientific_name) %>% 
    sample_n(1)
  #df %>% head %>% print
  # getting OTU
  otu = otu[,df$sample] %>% t 
  # subsampling tree
  to_rm = setdiff(tree$tip, rownames(otu))
  tree = drop.tip(tree, to_rm)
  # creating phylo4d 
  tree_4d = phylobase::phylo4d(tree, tip.data=otu)
  return(tree_4d)
}
#' lipamoran analysis on each trait
lipamoran_each = function(i, tree_4d, traits, reps=99){
  # phylosignal
  doParallel::registerDoParallel(threads)
  lipa_res = plyr::llply(as.list(traits), lipaMoran_per_OTU, tree_4d=tree_4d[[i]], reps=reps, .parallel=TRUE)
  # formatting
  lipa_res = do.call(rbind, lipa_res)
  lipa_res$subsample_rep = i
  return(lipa_res)
}

###read otu tabale and taxonomy
otu=read.csv("Relative abundance of gut microbial genera.csv",header = T,row.names = 1)
head(otu)
tax=read.csv("tax-genus.csv",header = T)
head(tax)

#filter
otu <- otu[which(rowSums(otu) >= 0.005), ]
host_tree = read.tree("Phylogenetic relationships of the host.nwk")
plot(host_tree)
any(host_tree$edge.length < 0)  # TRUE indicates that a negative value exists
which(host_tree$edge.length < 0)
host_tree$edge.length[host_tree$edge.length < 0]
host_tree$edge.length[host_tree$edge.length < 0] <- 0
host_tree$tip.label
# convert to phylo4d
host_tree_4d = phylobase::phylo4d(host_tree, tip.data=t(otu))
host_tree_4d %>% summary
ape::is.binary.tree(host_tree)
any(host_tree$edge.length == 0)
host_tree$edge.length[host_tree$edge.length == 0] <- 1e-6
# phylosignal calculation
physig_res = phylosignal::phyloSignal(host_tree_4d, methods = 'all', reps = 999)
# formatting output
tmp1 = physig_res$stat 
tmp1$OTU = rownames(tmp1)

tmp2 = physig_res$pvalue
tmp3 = apply(tmp2, 2, function(x) p.adjust(x, method='BH')) %>% as.data.frame
tmp2$OTU = rownames(tmp2)
tmp3$OTU = rownames(tmp2)

tmp1 = tmp1 %>%
  gather(method, coef, -OTU)
tmp2 = tmp2 %>%
  gather(method, pvalue, -OTU)
tmp3 = tmp3 %>%
  gather(method, qvalue, -OTU)

physeq_res_j = tmp1 %>%
  inner_join(tmp2, c('OTU', 'method')) %>%
  inner_join(tmp3, c('OTU', 'method'))

tmp1 %>% nrow %>% print
physeq_res_j %>% status
# plotting coef distribution
ggplot(physeq_res_j, aes(method, coef)) +
  geom_boxplot() +
  theme_bw()
# plotting pvalue distribution
ggplot(physeq_res_j, aes(coef, pvalue, color=method)) +
  geom_point(alpha=0.5) +
  theme_bw()
# plotting qvalue distribution
ggplot(physeq_res_j, aes(coef, qvalue, color=method)) +
  geom_point(alpha=0.5) +
  theme_bw()
# add taxonomy
physeq_res_j = physeq_res_j %>%
  inner_join(tax, c('OTU')) 
## all significant OTU
physeq_res_j_f = physeq_res_j %>%
  filter(pvalue < 0.05) %>%
  group_by(OTU) %>%
  mutate(n_methods = method %>% unique %>% length) %>%
  ungroup() %>%
  mutate(Phylum = Phylum %>% as.character,
         Phylum = Phylum %>% reorder(Domain %>% as.factor %>% as.numeric))
write.csv(physeq_res_j_f,"physeq_res_j_f.csv")
# plotting coef. values for each taxonomic group
ggplot(physeq_res_j_f, aes(Domain, coef, color=Domain)) +
  geom_boxplot() +
  facet_grid(method ~ ., scales='free_y') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=8)
  )

#####Phylosignal sensitivity analysis subsampling 1 sample per species#######
####local phylogenetic signal
lipaMoran_per_OTU = function(trait, tree_4d, reps=9999){
  res = phylosignal::lipaMoran(host_tree_4d, trait=traits, reps=9999, prox.phylo = "nNodes")
  z = colnames(res$lipa)[1]
  x = res$lipa
  colnames(x) = c('coef')
  y = res$p.value 
  colnames(y) = c('pvalue')
  df = cbind(x,y)
  df=as.data.frame(df)
  df$OTU = z
  df$host = rownames(df)
  rownames(df) = 1:nrow(df)
  return(df)
}
traits = physeq_res_j_f$OTU %>% unique
traits %>% length
# running in parallel
registerDoParallel(cores=2)
lipa_res = plyr::llply(as.list(traits),lipaMoran_per_OTU,tree_4d=host_tree_4d, reps=9999, .parallel=TRUE)
lipa_res = do.call(rbind, lipa_res)
# plotting 
ggplot(lipa_res, aes(coef, pvalue)) +
  geom_point(alpha=0.5) +
  theme_bw() 
# adjusting p-values
lipa_res$qvalue = p.adjust(lipa_res$pvalue, method='BH')
lipa_res$qvalue %>% summary
# significant Genera
lipa_res_f = lipa_res %>%
  filter(pvalue < 0.05) %>%
  inner_join(tax, c('OTU'))
# number of significant Genera
lipa_res_f$OTU %>% unique %>% length
write.csv(lipa_res_f,"lipa_res_f.csv")