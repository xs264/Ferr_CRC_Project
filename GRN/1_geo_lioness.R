library(tidyverse)

## ferroptosis genes
genes = read.csv("data/Ferrgene_GSE.csv")
genes = genes$Ferr_gene

## all lioness networks
all_files = list.files(path = "data/GEO_colon/", pattern = "*.csv", full.names = TRUE)

gene_score_all = matrix(NA,nrow = 12817,ncol = length(all_files))
TF_score_all = matrix(NA,nrow = 661,ncol = length(all_files))
gene_score_ft = matrix(NA,nrow = 344,ncol = length(all_files))
TF_score_ft = matrix(NA,nrow = 661,ncol = length(all_files))

for (i in 1:length(all_files)){
  
  print(i)
  file = all_files[i]
  
  ## full net node degrees
  net = fread(file)
  net = column_to_rownames(net,'TF') #661 * 12817
  TF.deg = rowSums(net)
  TG.deg = colSums(net)
  gene_score_all[,i] = TG.deg
  TF_score_all[,i] = TF.deg
  
  ## ferrop net node degrees
  net2 = net[,colnames(net)%in%genes]
  TF.deg = rowSums(net2)
  TG.deg = colSums(net2)
  gene_score_ft[,i] = TG.deg
  TF_score_ft[,i] = TF.deg
  
}

## add names
sample_names = sub(".*lioness_(.*)\\.csv", "\\1", all_files)

rownames(gene_score_all) = colnames(net)
colnames(gene_score_all) = sample_names

rownames(TF_score_all) = rownames(net)
colnames(TF_score_all) = sample_names

rownames(gene_score_ft) = intersect(colnames(net),genes)
colnames(gene_score_ft) = sample_names

rownames(TF_score_ft) = rownames(net)
colnames(TF_score_ft) = sample_names

## save rds
res = list("TG_score_all"=gene_score_all,
           "TF_score_all"=TF_score_all,
           "TG_score_ft"=gene_score_ft,
           "TF_score_ft"=TF_score_ft)


saveRDS(res,"result/geo_node_score.rds")


