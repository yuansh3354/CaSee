# clean
rm(list = ls())
gc()

# packages
library(Seurat)
library(stringr)
library(magrittr)
library("scales")
library(RColorBrewer)
library(ggplot2)
library(UCell)

# myfunctions 
myhead = function(x,DT=FALSE){
  print(paste('Data dimension: ',dim(x)[1],dim(x)[2]))
  if(dim(x)[2]>5){print(x[1:5,1:5])}
  else(print(head(x)))
  if(DT){DT::datatable(x,filter='top')}
}
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
)
show_col(mycolors)
# options
options(stringsAsFactors = F)
options(as.is = T)
df = read.csv('~/Desktop/1.csv')
library(ggplot2)
ggplot(data = df, aes(x = Sample_Origin, fill = Cell_type)) +
  geom_bar(position = "fill") +scale_fill_manual(values = mycolors[10:20]) +
  scale_y_reverse() +
  coord_flip()
ggsave2("Fig2A_barplot.pdf", path = "../results", width = 20, height = 5, units = "cm")
