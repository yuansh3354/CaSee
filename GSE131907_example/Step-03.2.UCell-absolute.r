# Calculate gene signature score
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
sce =readRDS('GSE131907_epi_std.rds')
gene_list = read.csv(file_path)
predict = read.csv('merge_data_h5_file/merge_data_predict_Integration_.csv',row.names = 1)
file_path = '/media/yuansh/14THHD/BscModel-V4/supplementary materials/DEGs/'
dir(file_path)
file_path = paste0(file_path,'GSE131907_Cancer_absolute_scale_wilcoxon.csv')
DEGs = read.csv(file_path,row.names = 1)
DEGs = DEGs[which( (DEGs$pvals_adj<0.05) & (abs(DEGs$logfoldchanges)>2)),]
############# 
sce = AddMetaData(sce,predict)
sce = sce[,which(!is.na(sce$Cell_subtype))]
sce = sce[,which(sce$Cell_subtype != 'Undetermined')]
sce$cnv.label = ifelse(sce$Cell_subtype %in% c('Malignant cells','tS1','tS2','tS3'),'Cancer','Normal')

sce$absolute = 'Others'
sce@meta.data[which((sce$cnv.label == 'Cancer') & (sce$scale_predict == 'Cancer')),'absolute'] ='Cancer'
sce@meta.data[which((sce$cnv.label == 'Normal') & (sce$scale_predict == 'Normal')),'absolute'] ='Normal'

module_names = colnames(gene_list)
module = list()
for(module_name in module_names){
  module[paste0(module_name,'_absolute')] = gene_list[,module_name] %>% unique() %>% 
           intersect(DEGs[,1]) %>% as.character() %>% list()
}
sce <- AddModuleScore_UCell(sce, features = module)

ids = sce@meta.data[which(sce$absolute !='Others'),] %>% rownames()
sce.plot = subset(sce,cells = ids)

paste0(module_names,'_absolute', "_UCell")
features = c('EMT_absolute_UCell','Energy.metabolism_absolute_UCell')

p1 = RidgePlot(sce.plot, features=features[1],cols = c("#E64B35B2",'#698EC3')
               , group.by = "absolute") 
p2 = RidgePlot(sce.plot, features=features[2],cols = c("#E64B35B2",'#698EC3')
               , group.by = "absolute") 
p3 = DotPlot(sce.plot,features=features,cols = c("lightgrey", "#e67d6b")
             ,group.by = "absolute")
x = CombinePlots(plots = list(p1,p2,p3),ncol = 1)
x
topptx(x,file="Module_score_absolute.pptx",  width = 12,
       height = 9)
