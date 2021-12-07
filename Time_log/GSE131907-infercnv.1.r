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