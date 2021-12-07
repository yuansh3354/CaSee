# TSNE
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
library(Seurat)
library(stringr)
library(magrittr)
library("scales")
library(RColorBrewer)
library(ggplot2)

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
sce = readRDS('rds/epi_raw.RDS')
predict = read.csv('merge_data_h5_file/merge_data_predict_Integration_.csv',row.names = 1)
identical(colnames(sce),rownames(predict))
sce = AddMetaData(sce,predict)

if(T){
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

nCount_lower <- 1000
nCount_upper <- 40000
nFeature_lower <- 500
nFeature_upper <- 7500
percent.mt_upper <- 18
percent.ncRNA_upper <-0.8
percent.LOC_upper <- 0.3

sce <- subset(sce, subset = nFeature_RNA > nFeature_lower & 
                    nFeature_RNA < nFeature_upper & 
                    nCount_RNA > nCount_lower & 
                    nCount_RNA < nCount_upper & 
                    percent.mt < percent.mt_upper & 
                    percent.ncRNA < percent.ncRNA_upper & 
                    percent.LOC < percent.LOC_upper
                    )

sce <- SCTransform(sce, verbose = T, vars.to.regress = c("nCount_RNA", "percent.mt",'percent.ncRNA','percent.LOC'), conserve.memory = T)
sce = RunPCA(object = sce, do.print = FALSE)

sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce)
sce <- RunUMAP(sce, dims = 1:4, verbose = T)
p1 = DimPlot(sce[,c(1:10)], reduction = "umap",group.by ='scale_predict',cols = c("#E64B35B2",'#698EC3')) 
p2 = DimPlot(sce[,c(1:10)], reduction = "umap",group.by ='Sample') + theme(legend.position="none")
x=CombinePlots(plots = list(p1,p2))
eoffice::topptx(x,'tsne_of_NC_cells.pptx',width = 10,height = 4)

p1 = DimPlot(sce, reduction = "umap",group.by ='scale_predict',cols = c("#E64B35B2",'#698EC3')) 
p2 = DimPlot(sce, reduction = "umap",group.by ='Sample') + theme(legend.position="none")
x=CombinePlots(plots = list(p1,p2))
ggsave("tsne_of_NC_cells.pdf", width = 20, height = 8, units = "cm")
