# clear
rm(list= ls())
library(magrittr)
setwd('/media/yuansh/14THHD/BscModel-V4/PMID33958794')
ids = dir(pattern  = 'cohort1')
ids

meta = read.csv(ids[1],row.names = 1)
expr = read.csv(ids[2],row.names = 1)

colnames(meta)
unique(meta$cellType)
colnames(expr)=gsub('\\.','-',colnames(expr)) 
use.cells = meta[which(meta$cell == 'Cancer_cell'),] %>% rownames()
cancer_count = expr[,use.cells]
write.csv(cancer_count,'cohort1_cancer_count.csv')


# clear
rm(list= ls())
library(magrittr)
setwd('/media/yuansh/14THHD/BscModel-V4/PMID33958794')
ids = dir(pattern  = 'cohort2')
ids

meta = read.csv(ids[1],row.names = 1)
expr = read.csv(ids[2],row.names = 1)

colnames(meta)
unique(meta$cellType)
colnames(expr)=gsub('\\.','-',colnames(expr)) 
use.cells = meta[which(meta$cell == 'Cancer_cell'),] %>% rownames()
cancer_count = expr[,use.cells]
write.csv(cancer_count,'cohort2_cancer_count.csv')