# Calculate gene signature score
### ---------------
###
### Create: Yuan.Sh
### Date: 2021-11-22 16:35:19
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Fujian Medical University
### Ref : https://github.com/carmonalab/UCell_demo
###
### ---------------

############# Step.00 ############# 
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
mycolors = c("#1a476f","#90353b","#55752f","#e37e00","#6e8e84",
             "#c10534","#938dd2","#cac27e","#a0522d","#7b92a8",
             "#2d6d66","#9c8847","#bfa19c","#ffd200","#d9e6eb",
             "#4DBBD5B2", "#00A087B2", "#3C5488B2", 
             "#F39B7FB2", "#8491B4B2","#91D1C2B2","#DC0000B2",
             "#7E6148B2","#E64B35B2",'#698EC3')
show_col(mycolors)
# options
options(stringsAsFactors = F)
options(as.is = T)
setwd('/media/yuansh/14THHD/BscModel-V4/GSE131907')

# main
expr = read.csv('GSE131907_epi_count.csv',row.names = 1)
meta = read.csv('GSE131907_epi_meta.csv',row.names = 1)
############# 
myhead(expr)
sce = CreateSeuratObject(counts = expr,meta.data = meta)
if(T){
  # Standerdize processing 
  ### 1.log
  sce = NormalizeData(object = sce,normalization.method =  "LogNormalize")
  ### 2.FindVariable
  sce = FindVariableFeatures(object = sce,selection.method = "vst", nfeatures = 2000)
  ### 3.ScaleData
  sce = ScaleData(object = sce)
  ### 4. PCA
  sce = RunPCA(object = sce, do.print = FALSE)
  ### 5.Neighbor
  sce= FindNeighbors(sce, dims = 1:20)
  ### 6. Clusters
  sce = FindClusters(sce) 
  ### 7.tsne
  sce=RunTSNE(sce,dims.use = 1:20)
  sce=RunUMAP(sce,dims = 1:20)
}
saveRDS(sce,file = 'GSE131907_epi_std.rds')
