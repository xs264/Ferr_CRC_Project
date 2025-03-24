library(data.table)
library(tidyverse)
tcga_lioness = fread("data/TCGA_colon/Colon_cancer_TCGA_AllSamples.csv")

## common samples
tcga = read.csv("data/TCGA-COADREAD.csv")
tcga_id = tcga$PATIENT_ID
lioness_id = colnames(tcga_lioness)[-(1:4)]
lioness_id = substr(lioness_id,1,12)
common_id = intersect(tcga_id,lioness_id) #346

## filter
tcga_net = tcga_lioness[,c(1,2,(which(lioness_id %in% common_id)+4)),with=F]
rm("tcga_lioness")

## Replace original sample IDs with common IDs
colnames(tcga_net)[-(1:2)] = substr(colnames(tcga_net)[-(1:2)],1,12)
saveRDS(as.data.frame(tcga_net),"data/tcga_net.rds")

##############################################################################
## calculate gene degrees
tcga_net = readRDS("data/tcga_net.rds")
tcga_net = as.data.table(tcga_net)
dt_Gene = tcga_net[, lapply(.SD, sum), by = gene, .SDcols = -c("TF", "gene")]
dt_TF = tcga_net[, lapply(.SD, sum), by = TF, .SDcols = -c("TF", "gene")]
gene_score_all = column_to_rownames(dt_Gene,"gene")
TF_score_all = column_to_rownames(dt_TF,'TF')

## ferr genes
genes = read.csv("data/Ferrgene_TCGA.csv")
genes = genes$Ferr_gene
tcga_net = tcga_net[tcga_net$gene %in% genes,]
tcga_net = as.data.table(tcga_net)

dt_Gene = tcga_net[, lapply(.SD, sum), by = gene, .SDcols = -c("TF", "gene")]
dt_TF = tcga_net[, lapply(.SD, sum), by = TF, .SDcols = -c("TF", "gene")]
gene_score_ft = column_to_rownames(dt_Gene,"gene")
TF_score_ft = column_to_rownames(dt_TF,'TF')

## save rds
res = list("TG_score_all"=gene_score_all,
           "TF_score_all"=TF_score_all,
           "TG_score_ft"=gene_score_ft,
           "TF_score_ft"=TF_score_ft)

saveRDS(res,"result/tcga_node_score.rds")
