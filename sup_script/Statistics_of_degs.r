###############################################################################################
setwd('/media/yuansh/14THHD/胶囊单细胞/Capsule/IPA/')
library(magrittr)
rm(list = ls())
dir_gene_statistic = NULL
ps = c(0.05,0.01,0.005,0.001)
fcs = c(0,1,1.5,2,2.5,3,3.5,4)
datasets = dir(pattern = '.csv')
for(i in 1:length(datasets)){
  assign(datasets[i], read.csv(datasets[i],header = T,row.names = 1))
}
for(i in 1:length(datasets)){
  ids = get(datasets[i])
  n_gene = NULL
  name = NULL
  
  for(p in ps){
    for (fc in fcs) {
      ids1 = which(ids$pvals_adj<p & abs(ids$logfoldchanges)>=fc) %>% length()
      n_gene = c(n_gene,ids1)
      name = c(name,paste0("p=",p,' & ','fc>',fc))
    }
  }
  ids2 = data.frame(n_gene,row.names = name) %>% set_colnames(datasets[i]) %>% as.matrix()
  dir_gene_statistic = cbind(dir_gene_statistic,ids2)
}
write.csv(dir_gene_statistic,'/media/yuansh/14THHD/胶囊单细胞/PRJNA591860_dif-statistic.csv')
###############################################################################################
setwd('/media/yuansh/14THHD/胶囊单细胞/测试集/GSE131907/IPA')
library(magrittr)
rm(list = ls())
dir_gene_statistic = NULL
ps = c(0.05,0.01,0.005,0.001)
fcs = c(0,1,1.5,2,2.5,3,3.5,4)
datasets = dir(pattern = '.csv')
for(i in 1:length(datasets)){
  assign(datasets[i], read.csv(datasets[i],header = T,row.names = 1))
}
for(i in 1:length(datasets)){
  ids = get(datasets[i])
  n_gene = NULL
  name = NULL
  
  for(p in ps){
    for (fc in fcs) {
      ids1 = which(ids$pvals_adj<p & abs(ids$logfoldchanges)>=fc) %>% length()
      n_gene = c(n_gene,ids1)
      name = c(name,paste0("p=",p,' & ','fc>',fc))
    }
  }
  ids2 = data.frame(n_gene,row.names = name) %>% set_colnames(datasets[i]) %>% as.matrix()
  dir_gene_statistic = cbind(dir_gene_statistic,ids2)
}
write.csv(dir_gene_statistic,'/media/yuansh/14THHD/胶囊单细胞/GSE131907_dif-statistic.csv')
###############################################################################################
setwd('/media/yuansh/14THHD/胶囊单细胞/GSE132465/IPA/')
library(magrittr)
rm(list = ls())
dir_gene_statistic = NULL
ps = c(0.05,0.01,0.005,0.001)
fcs = c(0,1,1.5,2,2.5,3,3.5,4)
datasets = dir(pattern = '.csv')
for(i in 1:length(datasets)){
  assign(datasets[i], read.csv(datasets[i],header = T,row.names = 1))
}
for(i in 1:length(datasets)){
  ids = get(datasets[i])
  n_gene = NULL
  name = NULL
  
  for(p in ps){
    for (fc in fcs) {
      ids1 = which(ids$pvals_adj<p & abs(ids$logfoldchanges)>=fc) %>% length()
      n_gene = c(n_gene,ids1)
      name = c(name,paste0("p=",p,' & ','fc>',fc))
    }
  }
  ids2 = data.frame(n_gene,row.names = name) %>% set_colnames(datasets[i]) %>% as.matrix()
  dir_gene_statistic = cbind(dir_gene_statistic,ids2)
}
write.csv(dir_gene_statistic,'/media/yuansh/14THHD/胶囊单细胞/GSE132465_dif-statistic.csv')