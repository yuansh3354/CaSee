# different gene expression counts in cancer and normal cells
### ---------------
###
### Create: Yuan.Sh
### Date: 2021-11-22 16:35:19
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Fujian Medical University
### Ref : DOI:https://doi.org/10.1016/j.cell.2020.07.017 (Figure S1F)
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
library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())

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
setwd('/media/yuansh/1THHD/Yuansh_Share/乳腺癌/data/data1/output')

############## Cohort1 ############## 
file_path = '/media/yuansh/My Passport/科研项目/深度学习模型/CapsuleNet预测肿瘤细胞/ref_data/Cancer_functional_module_gene_list.csv'
sce =readRDS('cohort1_std_analysis.rds')
gene_list = read.csv(file_path)
predict = read.csv('../cohort1_predict_Integration_.csv',row.names = 1)
file_path = '/media/yuansh/14THHD/BscModel-V4/PMID33958794/IPA-time/'
dir(file_path)
file_path = paste0(file_path,'scale_predict_wilcoxon_cohort1_Pre.csv')
DEGs = read.csv(file_path,row.names = 1)
DEGs = DEGs[which( (DEGs$pvals_adj<0.01) & (abs(DEGs$logfoldchanges)>1)),]

# using pre
sce = AddMetaData(sce,predict)
a = sce@meta.data
a = a[which(!is.na(a$scale_predict)),]
a$Log2_nCount_RNA = log(a$nCount_RNA,2)
df = a[which(a$timepoint == 'Pre'),]

p1 = ggplot(df,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('Cohort1_Pre')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")
df = a[which(a$timepoint == 'On'),]
p2 = ggplot(df,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('Cohort1_On')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which((a$timepoint == 'Pre') & (a$expansion == 'E')),]
p3 = ggplot(df,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot()+ ggtitle('Pre_E')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")
df = a[which((a$timepoint == 'Pre') & (a$expansion == 'NE')),]
p4 = ggplot(df,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot()+ ggtitle('Pre_NE')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which((a$timepoint == 'On') & (a$expansion == 'E')),]
p5 = ggplot(df,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot()+ ggtitle('On_E')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")
df = a[which((a$timepoint == 'On') & (a$expansion == 'NE')),]
p6 = ggplot(df,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot()+ ggtitle('On_NE')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")
x = CombinePlots(plots = list(p1,p2,p3,p4,p5,p6),legend='right')
topptx(x,filename = '~/Desktop/cohort1.pptx',width = 8,height = 6)

#wilcox.test
############## GSE131907 ############## 
setwd('/media/yuansh/14THHD/BscModel-V4/GSE131907/')
rm(list = ls())
gc()
predict = read.csv('merge_data_h5_file/merge_data_predict_Integration_.csv',row.names = 1)
sce = readRDS('rds/GSE131907_epi_std.rds')
sce = AddMetaData(sce,predict)

a = sce@meta.data
a$Log2_nCount_RNA = log(a$nCount_RNA,2)
p1 = ggplot(a,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('GSE131907')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which(a$Sample_Origin == 'mBrain'),]
p2 = ggplot(a,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('mBrain')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which(a$Sample_Origin == 'tLung'),]
p3 = ggplot(a,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('tLung')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which(a$Sample_Origin == 'tL/B'),]
p4 = ggplot(a,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('tL/B')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which(a$Sample_Origin == 'mLN'),]
p5 = ggplot(a,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('mLN')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

df = a[which(a$Sample_Origin == 'PE'),]
p6 = ggplot(a,aes(scale_predict,Log2_nCount_RNA,fill=scale_predict)) +geom_boxplot() + ggtitle('PE')+
  scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+ stat_compare_means(method = "t.test")

x = CombinePlots(plots = list(p1,p2,p3,p4,p5,p6),legend='right')
topptx(x,filename = '~/Desktop/GSE131907.pptx',width = 8,height = 6)
