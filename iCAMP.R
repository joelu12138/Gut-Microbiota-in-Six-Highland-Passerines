#loading package
library(iCAMP)
library(ape)
#
com.file="otus.csv"
tree.file="rep_phylo.tre"
treat.file="group.csv"
save.wd="C://Rworkplace"
clas<-"tax_table.csv"
#
rand.time=100
nworker=4
memory.G=50
#
setwd("C://Rworkplace")
comm=t(read.csv(com.file,header = T,row.names = 1))
tree=read.tree(file = tree.file)
treat=read.csv(treat.file,header = T,row.names = 1)
clas<-read.csv(clas,header = T,row.names = 1)
#
sampid.check=match.name(rn.list=list(comm=comm,treat=treat))
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
#
spid.check=match.name(cn.list=list(comm=comm),tree.list=list(tree=tree))
comm=spid.check$comm
tree=spid.check$tree
#The phylogenetic distance matrix was calculated
setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
}else{
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
  }

#iCAMP analysis
sig.index="Confidence"
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = 24, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
#
icbin=icamp.bins(icamp.detail = icres$detail,treat = treat,
                 clas=clas,silent=FALSE, boot = TRUE,
                 rand.time = rand.time,between.group = TRUE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda"))
write.csv(icbin$Pt,"ProcessImportance_EachGroup.csv")
write.csv(icbin$Ptk,"ProcessImportance_EachBin_EachGroup.csv")
write.csv(icbin$Ptuv,"ProcessImportance_EachTurnover.csv")
write.csv(icbin$BPtk,"BinContributeToProcess_EachGroup.csv")