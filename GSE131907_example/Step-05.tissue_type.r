# pie plot
### ---------------
###
### Create: Yuan.Sh
### Date: 2021-11-26 16:57:56
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Fujian Medical University
### Ref : 
###
### ---------------

############# Step.00 ############# 
# clean
rm(list = ls())
gc()

# packages
library("scales")
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
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
norm_range = function(x){
  ids = summary(x)
  iqr = 1.5*(ids[5] - ids[2])
  print(ids)
  print(paste0('norm_range:',ids[2]-iqr,'-',ids[5]+iqr))
}
show_col(mycolors)
# options
options(stringsAsFactors = F)
options(as.is = T)
setwd('/media/yuansh/14THHD/BscModel-V4/GSE131907')

# main
# expr = read.csv('GSE131907_epi_count.csv',row.names = 1)
# meta = read.csv('GSE131907_epi_meta.csv',row.names = 1)
# myhead(expr)
# sce = CreateSeuratObject(counts = expr,meta.data = meta)
# saveRDS(sce, file = "rds/epi_raw.RDS")
# remove(meta)
# remove(expr)
# gc()
# sce@meta.data$article_label = ifelse(sce$Cell_subtype %in% c('Malignant cells',"tS1",'tS2','tS3'),'Cancer','Normal')
# sce@meta.data$cnv.label = sce@meta.data$article_label
sce = readRDS('rds/epi_raw.RDS')
predict = read.csv('merge_data_h5_file/merge_data_predict_Integration_.csv',row.names = 1)
identical(colnames(sce),rownames(predict))
sce = AddMetaData(sce,predict)

if(F){
  erccs = grep('^ERCC-', x= rownames(sce),value = T) # value = T 获取名字
  rp = grep("^RP[SL][[:digit:]]", x= rownames(sce),value = T) # value = T 获取名字
  mt = grep('^MT-', x= rownames(sce),value = T) # value = T 获取名字) 
  # ncRNA = grep(".*\\.[0-9]", x= rownames(sce),value = T) # value = T 获取名字
  ncRNA = grep("^[A-Z][A-Z][0-9]*\\.[0-9]", x= rownames(sce),value = T)
  LOC = grep('(^LOC|LINC)[1-9]*', x= rownames(sce),value = T) # value = T 获取名字) 
  sce[["percent.ercc"]]  = PercentageFeatureSet(sce, pattern = "^ERCC-")
  sce[["percent.rp"]]  = PercentageFeatureSet(sce, pattern = "^RP[SL][[:digit:]]")
  sce[["percent.mt"]]  = PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ncRNA"]]  = PercentageFeatureSet(sce, pattern = ".*\\.[0-9]")
  sce[["percent.LOC"]]  = PercentageFeatureSet(sce, pattern = "(^LOC|LINC)[1-9]*")
  norm_range(sce@meta.data$nCount_RNA)
  norm_range(sce@meta.data$nFeature_RNA)
  norm_range(sce@meta.data$percent.mt)
  norm_range(sce@meta.data$percent.ncRNA)
  norm_range(sce@meta.data$percent.LOC)
}

df = table(sce@meta.data[,c('Sample_Origin','scale_predict')]) %>% as.data.frame()
theme_set(theme_cowplot())
x = ggplot(df, mapping = aes(Sample_Origin,Freq,fill = scale_predict)) +
  geom_bar(position="fill", stat = "identity")+
  theme(legend.title = element_blank())
topptx(x,'/media/yuansh/1THHD/Yuansh_Share/中转文件/2021年11月29日/sup_GSE131970_tiss.pptx')
table(sce@meta.data[,c('Sample_Origin','scale_predict')])
Cancer = table(sce@meta.data[,c('Sample_Origin','scale_predict')])[,1] / table(sce@meta.data$Sample_Origin)
Normal = table(sce@meta.data[,c('Sample_Origin','scale_predict')])[,2] / table(sce@meta.data$Sample_Origin)
write.csv(rbind(Cancer,Normal) %>% t,'/media/yuansh/14THHD/BscModel-V4/supplementary materials/旧版/GSE131970_tiss.csv')


sce@meta.data[which(is.na(sce$Cell_subtype)),'article_label'] = 'Undetermined'
sce@meta.data[which(sce$Cell_subtype == 'Undetermined'),'article_label'] = 'Undetermined'

df = table(sce@meta.data[,c('Sample_Origin','article_label')]) %>% as.data.frame()
theme_set(theme_cowplot())
x = ggplot(df, mapping = aes(Sample_Origin,Freq,fill = article_label)) +
  geom_bar(position="fill", stat = "identity")+scale_fill_manual(values = c("#E64B35B2",'#698EC3',"#7E6148B2"))+
  theme(legend.title = element_blank())
x
topptx(x,'/media/yuansh/1THHD/Yuansh_Share/中转文件/2021年11月29日/sup_GSE131970_article.pptx')

df = table(sce@meta.data[,c('Cell_subtype','scale_predict')]) %>% as.data.frame()
theme_set(theme_cowplot())
x = ggplot(df, mapping = aes(Cell_subtype,Freq,fill = scale_predict)) +
  geom_bar(position="fill", stat = "identity")
  theme(legend.title = element_blank())
x
topptx(x,'/media/yuansh/1THHD/Yuansh_Share/中转文件/2021年11月29日/sup_GSE131970_Cellsub_predict.pptx')



sum(sce$article_label == sce$scale_predict) / dim(sce@meta.data)[1]






