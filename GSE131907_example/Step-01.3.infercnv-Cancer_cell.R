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
if(T){
  load("./rawdata/sceRaw.Rdata")
  a = sce@meta.data
  table(a$Cell_type) %>% names
  cell.use = c("Epithelial cells")
  use.cells = a[which(a$Cell_type %in% cell.use),] %>% rownames
  sce <-subset(sce, cells=use.cells) 
}
if(T){
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
  predict = read.csv('/media/yuansh/14THHD/BscModel-V4/GSE131907/merge_data_h5_file/merge_data_predict_Integration_.csv',row.names = 1)
  sce = AddMetaData(sce,predict)
  sce@meta.data
  dfcount = as.data.frame(sce@assays$RNA@counts)
  # 第二个文件样本信息矩阵
  groupinfo= data.frame(cellId = colnames(dfcount))
  identical(groupinfo[,1],sce@meta.data$cell_id)
  groupinfo$cellType = sce@meta.data$scale_predict
  
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
  expFile='/media/yuansh/14THHD/BscModel-V4/GSE131907/9709_34383_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='/media/yuansh/14THHD/BscModel-V4/GSE131907/9709_34383_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='/media/yuansh/14THHD/BscModel-V4/GSE131907/9709_34383_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
}
names(table(groupinfo$cellType))
if(T){
  rm(list=ls())
  setwd('/media/yuansh/14THHD/BscModel-V4/GSE131907/')
  options(stringsAsFactors = F)
  library(Seurat)
  library(ggplot2)
  library(infercnv)
  expFile='9709_34383_expFile.txt' 
  groupFiles='9709_34383_groupFiles.txt'  
  geneFile='9709_34383_geneFile.txt'
  library(infercnv)
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c("Normal")) # 如果有正常细胞的话，把正常细胞的分组填进去
  #future::plan("multiprocess",workers=14)# 多核并行处理
  infercnv_all = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= "9709_34383_inferCNV",  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               num_threads=32,
                               denoise=F,
                               HMM=F)
}

library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "better_plot",output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色
