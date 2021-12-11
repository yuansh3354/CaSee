# Calculate the number of genes 
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
library(eoffice)
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
file_path = '/media/yuansh/My Passport/科研项目/深度学习模型/CapsuleNet预测肿瘤细胞/ref_data/Cancer_functional_module_gene_list.csv'
sce =readRDS('rds/GSE131907_epi_std.rds')
gene_list = read.csv(file_path)
predict = read.csv('merge_data_h5_file/merge_data_predict_Integration_.csv',row.names = 1)
file_path = '/media/yuansh/14THHD/BscModel-V4/supplementary materials/DEGs/'

############# 
sce = AddMetaData(sce,predict)
a = sce@meta.data

ggplot(a,aes(scale_predict,nFeature_RNA)) +geom_boxplot()


