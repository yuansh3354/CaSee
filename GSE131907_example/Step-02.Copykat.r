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
setwd('/media/yuansh/14THHD/胶囊单细胞/测试集/GSE131907/')

# main 
# 导入原始表达矩阵
# main 
# 导入原始表达矩阵
load("./rawdata/sceRaw.Rdata")

a = sce@meta.data
table(a$Cell_type) %>% names
cell.use = c("Endothelial cells","Epithelial cells","Fibroblasts")
use.cells = a[which(a$Cell_type %in% cell.use),] %>% rownames
sce <-subset(sce, cells=use.cells) 
library(copykat)
if(F){
  # because this matrix is too large to run, so we decided to rm some genes
  RP = grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = sce), value = FALSE) # 获取 index
  sce = sce[-RP,]
  MT = grep(pattern = "^MT-", x = rownames(x = sce), value = FALSE) # 获取 index
  sce = sce[-MT,]
  
  
  sce <- NormalizeData(object = sce, scale.factor = 1e6)
  sce <- FindVariableFeatures(object = sce,nfeatures=15000)
  
  ids = sce@assays$RNA@var.features
  sce = sce[ids,]
  
}
start = Sys.time()
print(start)
exp.rawdata = sce@assays$RNA@counts
copykat.test <- copykat(rawmat=exp.rawdata, 
                        id.type="S", 
                        cell.line="no", 
                        ngene.chr=15, 
                        win.size=25, 
                        KS.cut=0.15, 
                        sam.name="GSE131907", 
                        distance="euclidean", 
                        n.cores=16)
finish= Sys.time()
pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

