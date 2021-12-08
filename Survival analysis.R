# Survival analysis
### ---------------
###
### Create: Yuan.Sh
### Date: 2021-11-24 16:26:24
### Email: yuansh3354@163.com
### Blog: https://blog.csdn.net/qq_40966210
### Fujian Medical University
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
library(scales)
library(RColorBrewer)
library(ggplot2)
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
setwd('/media/yuansh/14THHD/BscModel-V4/supplementary materials')

############# ANKRD30A.tsv ############# 
df1 = read.table('ANKRD30A.tsv',sep = '\t',header = T)
df1$gene.expression.RNAseq...HTSeq...FPKM..ANKRD30A = ifelse(df1$data < 0.6 ,'low','high')
library(survival)
library(survminer)
library(magrittr)
df = df1[,c(1,4,3,5)]

if(T){
  annotate_x <- 20
  plot_title_name = NULL
  #df = cbind(rawdata, label)
  survival_table=df
  meta = survival_table
  colnames(meta) = c('ID','time',"event",'group')
  meta = na.omit(meta)
  sfit1 <- survfit(Surv(time, event)~group, data=meta) #primary_ER/CTC_ER
  summary = summary(coxph(Surv(time, event)~group, data=meta))
  HR = round(summary$coefficients[2],2)
  CI95 = round(summary$conf.int[3:4],2)
  LogRankP = signif(summary$sctest[3],digits = 3)
  pp = paste("LogRank p = ",LogRankP)
  HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
  x= ggsurvplot(
    sfit1,palette =c("#E64B35B2",'#698EC3'),
    pval =F,
    font.y=14,
    pval.size = 8,
    pval.coord = c(0.2,0.1),
    conf.int = F,
    risk.table ='nrisk_cumcensor',
    #legend.labs =c("Low Rate", "High Rate"), 
    xlab ="Time[months]", 
    ylab ="Overall Survival",
    #legend.title = "PDL-1 Positive+++ rate", 
    legend = c(0.9,0.9),
    risk.table.title = 'Number at Risk \n(number censored)'
  )+theme_survminer(
    base_family = "Times New Roman",
    base_size = 14,
    font.main = c(14, "bold",'black'),
    font.submain = c(14, "bold",'black'),
    font.caption = c(14, "bold",'black'),
    font.x = c(14, "bold",'black'),
    font.y = c(14, "bold",'black'),
    font.tickslab = c(14, "bold",'black'),
    font.legend = c(14, "bold",'black')
  )
}
x = CombinePlots(plots = list(x$plot,x$table),ncol = 1)
topptx(x,file="/media/yuansh/1THHD/Yuansh_Share/中转文件/2021年11月24日/ANKRD30A.pptx")

############# SYTL2.tsv ############# 

df1 = read.table('SYTL2.tsv',sep = '\t',header = T)
summary(df1$data)
#df1$gene.expression.RNAseq...HTSeq...FPKM..SYTL2 = ifelse(df1$data < 3.5 ,'low','high')
library(survival)
library(survminer)
library(magrittr)
df = df1[,c(1,4,3,5)]

if(T){
  annotate_x <- 20
  plot_title_name = NULL
  #df = cbind(rawdata, label)
  survival_table=df
  meta = survival_table
  colnames(meta) = c('ID','time',"event",'group')
  meta = na.omit(meta)
  sfit1 <- survfit(Surv(time, event)~group, data=meta) #primary_ER/CTC_ER
  summary = summary(coxph(Surv(time, event)~group, data=meta))
  HR = round(summary$coefficients[2],2)
  CI95 = round(summary$conf.int[3:4],2)
  LogRankP = signif(summary$sctest[3],digits = 3)
  pp = paste("LogRank p = ",LogRankP)
  HHRR = paste("HR = ",HR,"( 95% CI:",CI95[1],"-",CI95[2],")")
  x=ggsurvplot(
    sfit1,palette =c("#E64B35B2",'#698EC3'),
    pval =F,
    font.y=14,
    pval.size = 8,
    #pval.coord = c(0.2,0.1),
    conf.int = F,
    risk.table ='nrisk_cumcensor',
    #legend.labs =c("Low Rate", "High Rate"), 
    xlab ="Time[months]", 
    ylab ="Overall Survival",
    #legend.title = "PDL-1 Positive+++ rate", 
    legend = c(0.9,0.9),
    risk.table.title = 'Number at Risk \n(number censored)'
  )+theme_survminer(
    base_family = "Times New Roman",
    base_size = 14,
    font.main = c(14, "bold",'black'),
    font.submain = c(14, "bold",'black'),
    font.caption = c(14, "bold",'black'),
    font.x = c(14, "bold",'black'),
    font.y = c(14, "bold",'black'),
    font.tickslab = c(14, "bold",'black'),
    font.legend = c(14, "bold",'black')
  )
}
x = CombinePlots(plots = list(x$plot,x$table),ncol = 1)
topptx(x,file="/media/yuansh/1THHD/Yuansh_Share/中转文件/2021年11月24日/SYTL2.pptx")
