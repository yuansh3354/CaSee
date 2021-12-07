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
setwd('/media/yuansh/14THHD/胶囊单细胞/GSE132465/')

# main 
# 导入原始表达矩阵
load("processfile/GSE132465_seurat.RData")
sce = main_tiss_filtered

library(copykat)
start = Sys.time()
print(start)
exp.rawdata = sce@assays$RNA@counts
copykat.test <- copykat(rawmat=exp.rawdata, 
                        id.type="S", 
                        cell.line="no", 
                        ngene.chr=5, 
                        win.size=25, 
                        KS.cut=0.15, 
                        sam.name="TNBC1", 
                        distance="euclidean", 
                        n.cores=16)
finish= Sys.time()
pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

