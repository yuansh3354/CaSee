if(T){
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
  # 标准流程
  if(T){
    sce <- NormalizeData(object = sce, scale.factor = 1e6)
    sce <- FindVariableFeatures(object = sce)
    sce <- ScaleData(object = sce)
    sce <- RunPCA(object = sce, do.print = FALSE)
    n.pcs <- 20
    sce <- FindNeighbors(object = sce, dims = 1:20, verbose = T)
    sce <- FindClusters(object = sce, verbose = T, resolution = 0.5)
    sce <- RunTSNE(sce, dims = 1:20)
    p = DimPlot(sce, reduction = "tsne", label = TRUE)
    p
  }
}
epi = c(2,9,11,12,16,18,19,24)
cells.use.epi <- row.names(sce@meta.data)[which(sce@meta.data$RNA_snn_res.0.5 %in% epi)]
sce@meta.data$cell_annotion = 'Normal'
sce@meta.data[which(row.names(sce@meta.data)%in%cells.use.epi),'cell_annotion'] = 'epi'
table(sce@meta.data$cell_annotion)
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
  # 输出
  expFile='./processfile/23479_samples_expFile.txt'
  write.table(dfcount ,file = expFile,sep = '\t',quote = F)
  groupFiles='./processfile/23479_samples_groupFiles.txt'
  write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
  geneFile='./processfile/23479_samples_geneFile.txt'
  write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
}
names(table(groupinfo$cellType))
if(T){
  rm(list=ls())
  options(stringsAsFactors = F)
  library(Seurat)
  library(ggplot2)
  library(infercnv)
  expFile='./processfile/23479_samples_expFile.txt' 
  groupFiles='./processfile/23479_samples_groupFiles.txt'  
  geneFile='./processfile/23479_samples_geneFile.txt'
  library(infercnv)
  gc()
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                      annotations_file=groupFiles,
                                      delim="\t",
                                      gene_order_file= geneFile,
                                      ref_group_names=c('Normal')) # 如果有正常细胞的话，把正常细胞的分组填进去
  # Run infer CNV
  infercnv_all = infercnv::run(infercnv_obj,
                               cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= "./processfile/23479_samples_inferCNV",  # dir is auto-created for storing outputs
                               cluster_by_groups=F,   # cluster
                               num_threads=32,
                               denoise=F,
                               HMM=F)
}




