# clear
rm(list= ls())
library(magrittr)
setwd('~/Desktop/')
floders = dir(pattern = 'Normal_')
# importing Datasets


floder = floders[3]
floder
list = dir(floder,pattern = "Anno")
info = gsub('.csv','',list)
# using assing function to import dataset one-time
for(i in 1:length(info)){
  assign(info[i], read.csv(paste0(floder,'/',list[i]),sep = ',',header = T,row.names = 1))
}
list = dir(floder,pattern = "dge")
ids = gsub('.txt','',list)
ids = grep('dge',ids,value = T)
# using assing function to import dataset one-time
for(i in 1:length(ids)){
  print(list[i])
  assign(ids[i], read.csv(paste0(floder,'/',list[i]),sep = ',',header = T,row.names = 1))
}


i=5
samples = grep('_Anno',ls(),value = T) %>% gsub('_Anno','',.)
length(samples)
sample = samples[i]
count = grep(sample,ls(),value = T) %>% grep('_dge',.,value = T) %>% get()
info = grep(sample,ls(),value = T) %>% grep('_Anno',.,value = T) %>% get()
unique(info[,3])
epi = grep('Parietal|Gastric|Pit|Epithelial|AT2|AT1|Clara|epithelial|Hepatocyte|hepatocyte',info[,dim(info)[2]],value = T)
unique(epi)

info = info[which(info[,dim(info)[2]]%in% epi),] %>% rownames()
length(info)
count = count[,info]
dim(count)
sample
paste0(floder,'/',sub('.txt','.csv',list[i]))
write.csv(count,paste0(floder,'/',sub('.txt','.csv',list[i])))

# Proliferating
my_rm = grep('Lung',ls(),value = T)
rm(list = my_rm)
