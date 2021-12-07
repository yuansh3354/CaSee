# Running inferCNV
### ---------------
###
### Create: Yuan.Sh
### Date: 2021-12-01 11:07:17
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Fujian Medical University
### Ref : https://blog.csdn.net/qq_40966210/article/details/120543312?ops_request_misc=%257B%2522request%255Fid%2522%253A%2522163832806116780265428782%2522%252C%2522scm%2522%253A%252220140713.130102334.pc%255Fblog.%2522%257D&request_id=163832806116780265428782&biz_id=0&utm_medium=distribute.pc_search_result.none-task-blog-2~blog~first_rank_v2~rank_v29-2-120543312.pc_v2_rank_blog_default&utm_term=infer&spm=1018.2226.3001.4450
###
### ---------------

if(T){
  rm(list=ls())
  setwd('/media/yuansh/14THHD/BscModel-V4/GSE131907/script')
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
                               num_threads=32,
                               denoise=F,
                               HMM=F)
}