# clear
rm(list=ls())
# options
options(stringsAsFactors = F)
options(as.is = T)
# packages
library(stringr)
library(magrittr)
library(ggplot2)
library(h5)
library("org.Hs.eg.db")

# myfunctions 
myhead = function(x){
  if(dim(x)[2]>5){return(x[1:5,1:5])}
  else(head(x))
}
myscale = function(x, scale_factor){
  (x / sum(x)) *scale_factor
}

# set work-dir
# pass

# load data
df = read.table('TcgaTargetGtex_gene_expected_count', head=T, row.names=1)
tcga = read.csv('pan_cancer_annotion.csv', head=T,  row.names=1)
gtex = read.csv('gtex_annotion.csv', head=T,  row.names=1)

#save(df,tcga, gtex, file='data.rdata')
load('data.rdata')
# data-preprocess
myhead(df)
colnames(df) = gsub('\\.', '-',colnames(df))

# get samples 
sample.list = colnames(df)


# filter gtex samples
sample.group=str_split(colnames(df),'-',simplify = T)[,1] 
sample.list = sample.list[grep('GTEX|TCGA',sample.group)]
df = df[,sample.list]

# change sample ids
t = str_split(sample.list,'-',simplify = T) 
sample.list = paste(t[,1],t[,2],t[,3],t[,4],sep='-')

colnames(df) = sample.list
tcga = tcga[!duplicated(tcga$sampleID),]
gtex = gtex[!duplicated(gtex$SAMPID),]
rownames(tcga) = tcga$sampleID
rownames(gtex) = gtex$SAMPID

# get tcga samples and gtex samples
t = grep('TCGA', sample.list, value = T)
g = grep('GTEX', sample.list, value = T)
t = intersect(t, rownames(tcga))
g = intersect(g, rownames(gtex))
sample.list = c(t,g)
tcga = tcga[t,]
gtex = gtex[g,]

# get info
gtex$label = 0
colnames(gtex) = colnames(tcga)
meta = rbind(gtex,tcga) %>% na.omit()
myhead(meta)

ids = intersect(rownames(meta),colnames(df))
meta = meta[ids,]
raw = df
df = df[,ids]
identical(rownames(meta), colnames(df))
colnames(meta) = c("sampleID","Tissue","label")

# get gene symbols
rownames(df) = gsub('\\.[[:digit:]]','',rownames(df))
genes = mapIds(org.Hs.eg.db,keys=rownames(df),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
genes = genes[which(!is.na(genes))]
ensg = names(genes)
genes = as.character(genes)
df = df[ensg,]

# remove duplicated genes ~30min
tab = data.frame(ensg,genes)
df <- t(sapply(split(tab[,1], tab[,2]), function(ids){
  colMeans(df[ids,,drop=FALSE])
}))

# main
# check data 
myhead(df)
myhead(meta)
dim(df)
dim(meta)
identical(rownames(meta), colnames(df))
which(duplicated(rownames(df)))
which(duplicated(colnames(df)))
expr = df
save(meta, expr, file='gtex-tcga.rdata')
rm(list = ls())
# reload
load('gtex-tcga.rdata')
myhead(expr)
myhead(meta)
dim(expr)
dim(meta)
identical(rownames(meta), colnames(expr))
which(duplicated(rownames(expr)))
which(duplicated(colnames(expr)))
write.table(expr, file = 'expr.tsv')
write.table(meta, file = 'meta.tsv')
# and the load in python to get h5 files
"
Best Regards,
Yuan.SH
---------------------------------------
School of Basic Medical Sciences,
Fujian Medical University,
Fuzhou, Fujian, China.
please contact with me via the following ways:
(a) e-mail :yuansh3354@163.com
"
