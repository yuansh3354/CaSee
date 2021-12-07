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
  library(magrittcr)
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
if(F){
  load("./rawdata/sceRaw.Rdata")
  
  a = sce@meta.data
  table(a$Cell_type) %>% names
  cell.use = c("Endothelial cells","Epithelial cells","Fibroblasts")
  use.cells = a[which(a$Cell_type %in% cell.use),] %>% rownames
  sce <-subset(sce, cells=use.cells) 
}
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
# infercnv流程
if(F){
  # 第一个文件count矩阵
  dfcount = as.data.frame(sce@assays$RNA@counts)
  # 第二个文件样本信息矩阵
  groupinfo= data.frame(cellId = colnames(dfcount))
  identical(groupinfo[,1],sce@meta.data$cell_id)
  groupinfo$cellType = sce@meta.data$Cell_type
  
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
  expFile='./processfile/40231_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='./processfile/40231_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='./processfile/40231_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
}
names(table(groupinfo$cellType))
if(T){
  rm(list=ls())
  setwd('/media/yuansh/14THHD/胶囊单细胞/测试集/GSE131907/')
  options(stringsAsFactors = F)
  library(Seurat)
  library(ggplot2)
  library(infercnv)
  expFile='./processfile/40231_expFile.txt' 
  groupFiles='./processfile/40231_groupFiles.txt'  
  geneFile='./processfile/40231_geneFile.txt'
  library(infercnv)
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c("Endothelial cells","Fibroblasts")) # 如果有正常细胞的话，把正常细胞的分组填进去
  #future::plan("multiprocess",workers=14)# 多核并行处理
  infercnv_all = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= "./processfile/40231_inferCNV",  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               num_threads=28,
                               denoise=F,
                               HMM=F)
}




