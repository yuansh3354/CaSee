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
rm(main_tiss_filtered)
sce
# infercnv流程
sce@meta.data$Cell_type %>% table() %>% names()
sce@meta.data$cell_annotion = ifelse(sce@meta.data$Cell_type == "Epithelial cells",'epi','Normal')
sce@meta.data$cell_annotion %>% table()
if(T){
  # 第一个文件count矩阵
  dfcount = as.data.frame(sce@assays$RNA@counts)
  # 第二个文件样本信息矩阵
  groupinfo= data.frame(cellId = colnames(dfcount))
  identical(groupinfo[,1],sce@meta.data$cell_id)
  groupinfo$cellType = sce@meta.data$cell_annotion
  
  # 第三文件
  library(AnnoProbe)
  geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
  geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
  geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
  
  ## 这里可以去除性染色体
  # 也可以把染色体排序方式改变
  dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
  dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
  
  myhead(dfcount)
  myhead(geneInfor)
  myhead(groupinfo)
  # 输出
  expFile='./processfile/40922_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='./processfile/40922_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='./processfile/40922_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
}
names(table(groupinfo$cellType))

if(T){
  setwd('/media/yuansh/14THHD/胶囊单细胞/GSE132465/')
  rm(list=ls())
  options(stringsAsFactors = F)
  library(Seurat)
  library(ggplot2)
  library(infercnv)
  expFile='./processfile/40922_expFile.txt' 
  groupFiles='./processfile/40922_groupFiles.txt'  
  geneFile='./processfile/40922_geneFile.txt'
  library(infercnv)
  gc()
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c('Normal')) # 如果有正常细胞的话，把正常细胞的分组填进去
  # Run infer CNV
  infercnv_all = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= "./processfile/40922_inferCNV",  # dir is auto-created for storing outputs
                               cluster_by_groups=TRUE,   # cluster
                               num_threads=32,
                               denoise=F,
                               HMM=F)
}
