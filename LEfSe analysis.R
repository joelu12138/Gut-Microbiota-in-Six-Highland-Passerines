#loading package
library(microeco) # Microbial Community Ecology Data Analysis
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(magrittr) # A Forward-Pipe Operator for R

feature_table <- read.csv('feature_table.csv', row.names = 1)
sample_table <- read.csv('sample_table.csv', row.names = 1)
tax_table <- read.csv('tax_table.csv', row.names = 1)

head(feature_table)[1:6,1:6]; head(sample_table)[1:6, ]; head(tax_table)[,1:6]

# Create a microtable object
dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax_table)
dataset

# performed lefse analysis
lefse <- trans_diff$new(dataset = dataset, 
                        method = "lefse", 
                        group = "Group", 
                        alpha = 0.05,
                        p_adjust_method = "none",
                        lefse_subgroup = NULL,
)
data2=lefse$res_diff
write.csv(data2,"data.csv")
# From v0.8.0, threshold is used for the LDA score selection.
lefse$plot_diff_bar(threshold = 4)
# clade_label_level 5 represent phylum level in this analysis
p<-lefse$plot_diff_cladogram(
  clade_label_level = 6, 
  only_select_show = FALSE, 
  group_order = c("GRTI", "ROSP", "BLRS", "TWIT", "WRSF", "ETSP"),
  color =c('#0c4e9b','#f58693','#c72228','#ffbc80','#6b98c4','#f98f34'),
  node_size_offset = 1.5,
  annotation_shape = 21,# 处理图例不会一起变，仍是矩形。
  annotation_shape_size = 4.5,
  clade_label_size =1,
  alpha = 0.2,
  branch_size=5
)
ggsave(p,filename = "PCoA.pdf",width = 12, height = 8)