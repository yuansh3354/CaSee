# CaSee <img src="yuansh-logo.jpg" width="280px" align="right" />

The most important things for single cell RNA sequencing of human cancer is to identify which is cancer or normal cell. However in the present, the common method is to evaluate copy number variation (CNV). These methods which are based on the principle that CNV is common in human cancers about 85%.  And these methods are based on expression data to inferred CNV, so there is a certain error (~20%) with actual genome copy number compilation. In conclusion, using CNV to define cancer cells is time-consuming and inaccurate. Therefore, we developed CaSee, a tumor cell discrimination model based on deep learning framework. it's important to note that we developed this model only for identification of the cancer cells of scRNA and connot replace the copy number model. CaSee has stronger performance and better professionalism in cancer cell idetification.

## Pre-requisites:

- Linux (Based on Ubuntu 20.04 LTS, Personal Computer) 
- CPU AMD Ryzen 9 3950X
- NVIDIA GeForce RTX 3090 24GB 384bit 1695MHz 19500MHz 
- Memory 128G (32GB*4) DDR4 3200MHz

### Environment and resource allocation

---

For instructions on installing anaconda on your machine (download the distribution that comes with python 3):
https://www.anaconda.com/distribution/

Next, use the environment configuration file located in **configs/CaSee_env_info.yaml** to create a conda environment:

**It will take several minutes or hours depending on the speed of the network. Please be patient.**

```
conda env create -n CaSee -f configs/CaSee_env_info.yaml
```

## Prepare model args

The args config file is **configs/CaSee_Model_configs.yaml**

```
# configs of trainning-loop, plz attention this config_file must be in floder named "configs".
# 训练集使用的配置文件，请务必放在configs文件夹下 
# Date: 2021-11-30 14:32:03
# Author: Yuan.Sh
--- 
# training model args
data_arguments: 
  
  # Your data should be stored in work_dir, it can be raw count matrix or candidate cancer cell count matrix
  # The format of matrix must be csv, and the row is genesymbol, the col is cell_id
  work_dir: /media/yuansh/14THHD/BscModel-V4/GSE131907/ 
  
  # Your input file name
  Counts_expr: GSE131907_epi_count.csv # must .csv
  
  # You can give the Tissue_type as you know, Don't worry that the model will use this infomation. this information just use in weight the Probability, if you only use scale Probability, it will not be use.
  Tissue_type: "Tumor" # Tissue type come from ['Tumor','Adjacent','Normal','Unknow']
  
  # when you do not know the marker genes of candidate cancer cells, you can use "elimination method" to get cells. You just set use_others=True
  use_others: False # if you want to use cell cluseter which is not be annotation.
  
  # If you just want to use mRNA in step of cell_annotion, you can set remove_genes=True 
  remove_genes: False # remove mt, rp, ercc, LNC, non-coding RNA

# If you input file is raw count matrix, cell_annotion=True, else cell_annotion=False
# Because CaSee is not applicable in non-candidate cells 
cell_annotion: False    

# This is default cells markers, you can use your own marker to define cell cluster.
Marker_genes:
  T_cell: ["CD3D",'CD3E','CD2']
  Fibroblast: ['COL1A1','DCN','C1R']
  Myeloid: ['LYZ','CD68','TYROBP']
  B_cell: ['CD79A','MZB1','MS4A1']
  Endothelial: ['CLDN5','FLT1','RAMP2']
  Mast: ['CPA3','TPSAB1','TPSB2']
  DC: ['LILRA4','CXCR3','IRF7']
  Cancer: ['EPCAM']
  
save_files:
  files: merge_data

# if traing times==1, the p-val will not give. 
# if your expr's gene much different from the referenc matrix. you can choose suitable args of batch_size or lr
trainig_loop:
  times: 10
  batch_size: 128  # if GPU not enough, plz set sutiable number  
  max_epochs: 20 # Generally speaking 20 epochs is enough to get the result, if can't plz set sutiable number
  seed: 42
  lr: 0.0005 
  split_data_seed: 0
  
ckpt: 

```

## Running Model

```
python CaSee.py --config configs/CaSee_Model_configs.yaml
```

## Contact

```
School of Basic Medical Sciences,
Fujian Medical University,
Fuzhou, Fujian, China.
please contact with me via the following ways:
(a) e-mail :yuansh3354@163.com

---------------------------------------
Best Regards,
Yuan.SH
---------------------------------------
```

