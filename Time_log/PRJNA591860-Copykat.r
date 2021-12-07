# clear
if(T){
  rm(list=ls())
}

# options
if(T){
  options(stringsAsFactors = F)
  options(as.is = T)
}

# packages
if(T){
  library(stringr)
  library(magrittr)
  library(ggplot2)
  library(Seurat)
  library(devtools)
  library(clustree)
  library(tidyverse)
  library(gridExtra)
  library(ggridges)
  library(ggplot2)
  library(ggExtra)
  
}

# myfunctions 
if(T){
  myhead = function(x){
    if(dim(x)[2]>5){return(x[1:5,1:5])}
    else(head(x))
  }
  myscale = function(x, scale_factor){
    (x / sum(x)) *scale_factor
  }
}

# set work-dir
setwd('/media/yuansh/14THHD/胶囊单细胞/Capsule')

# main 
# 导入原始表达矩阵
load("./processfile/S01-osi-Seurat.RData")
sce = merge.data
rm(merge.data)
if(T){
  sce@meta.data$sample_name <- as.character(sce@meta.data$sample_name)
  sample_name <- as.character(sce@meta.data$sample_name)
  # Make table 
  tab.1 <- table(sce@meta.data$sample_name) 
  # Which samples have less than 10 cells 
  samples.keep <- names(which(tab.1 > 10))
  metadata_keep <- filter(sce@meta.data, sample_name %in% samples.keep)
  # Subset Seurat object 
  sce <- subset(sce, cells=as.character(metadata_keep$cell_id))
  sce
  table(sce@meta.data$sample_name)
  table(sce@meta.data$patient_id)
  #save(sce, file=paste(dir,"S03_subset_preprocessed.RData", sep=""))
}# 剔除组织样本中少于10个细胞的组织
raw = sce
rm(sce)

# # 导入注释信息
# pre.label = read.csv("./processfile/predict.csv",row.names = 1)
# # 提取非免疫细胞
# use.cells  <- rownames(pre.label)
# cnv.cells = raw@meta.data[which(!is.na(raw@meta.data$cnv.label)),] %>% rownames()
# use.cells = intersect(use.cells,cnv.cells)
# sce <-subset(raw, cells=use.cells)  # 从原始数据中获取非免疫细胞表达矩阵
# pre.label = pre.label$label %>% as.data.frame() %>% set_rownames(rownames(pre.label))
# sce = AddMetaData(sce,metadata = pre.label,col.name = 'prelabel')
# # sce@meta.data$prelabel = ifelse(is.na(sce@meta.data$prelabel),sce@meta.data$immune_annotation,sce@meta.data$prelabel)
# # copyKAT
if(F){
  sce = raw
  library(copykat)
  exp.rawdata = sce@assays$RNA@counts
  start = Sys.time()
  print(start)
  copykat.test <- copykat(rawmat=exp.rawdata, 
                          id.type="S", 
                          cell.line="no", 
                          ngene.chr=5, 
                          win.size=25, 
                          KS.cut=0.15, 
                          sam.name="Capsule", 
                          distance="euclidean", 
                          n.cores=16)
  finish= Sys.time()
  pred.test <- data.frame(copykat.test$prediction)
  # CNA.test <- data.frame(copykat.test$CNAmat)
  sce = AddMetaData(sce,metadata = pred.test,col.name = 'pred.test')
  save(pred.test,sce, file='processfile/copykat.Rdata')
}
# 13.18
copykat = read.csv('copykat_prediction.txt',sep = '\t',row.names = 1)

sce = AddMetaData(sce,metadata = copykat,col.name = 'copykat')

# check C-index 
sce@meta.data$copykat %>% unique()
sce@meta.data$cnv.label %>% unique()
sce@meta.data$prelabel %>% unique()
ids1 = sce@meta.data[which(sce@meta.data$prelabel !='Unknow'),] %>% rownames()
ids2 = sce@meta.data[which(!is.na(sce@meta.data$copykat)),] %>% rownames()
use.cells = intersect(ids1,ids2)
sce <-subset(sce, cells=use.cells)  
# 
sce@meta.data$copykat %>% unique()
sce@meta.data$cnv.label %>% unique()
sce@meta.data$prelabel %>% unique()

#
temp = sce@meta.data[,c('copykat','cnv.label','prelabel')]
temp$copykat = ifelse(temp$copykat == "aneuploid","Tumor", "Normal")
temp$cnv.label = ifelse(temp$cnv.label == "cancer cell","Tumor", "Normal")

which(temp$copykat == temp$cnv.label) %>% length()
which(temp$copykat == temp$prelabel) %>% length()
which(temp$cnv.label == temp$prelabel) %>% length()

