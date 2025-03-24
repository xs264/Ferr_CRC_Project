library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
tcga_limma = function(data,tcga_pheno){
  ## clinical variable
  mut = tcga_pheno$KRAS_mut
  mut[mut=="0: WT"] = "0"
  mut[mut=="1: Mutant"] = "1"
  mut = factor(mut)
  
  sex = factor(tcga_pheno$SEX)
  age = as.numeric(tcga_pheno$AGE)
  stage = tcga_pheno$STAGE
  stage[stage=="STAGE I"] = "I"
  stage[stage=="STAGE II"] = "II"
  stage[stage=="STAGE III"] = "III"
  stage = as.factor(stage)
  
  ### 1. differential targeting
  design = model.matrix(~sex + mut + sex*mut + age + stage)
  #design = model.matrix(~sex + mut + sex*mut)
  colnames(design) <- make.names(colnames(design))
  colnames(design)
  fit = lmFit(as.matrix(data),design)
  contrast = makeContrasts(mut_vs_wt_in_F = mut1,
                           mut_vs_wt_in_M = mut1 + sexMale.mut1,
                           levels = colnames(design))
  fit = contrasts.fit(fit, contrast)
  DEG_fit = eBayes(fit)
  
  DEG1 = topTable(DEG_fit,number=Inf,coef = "mut_vs_wt_in_F",
                  sort.by = "p",p.value = 1,lfc = 0)
  DEG2 = topTable(DEG_fit,number=Inf,coef = "mut_vs_wt_in_M",
                  sort.by = "p",p.value = 1,lfc = 0)
  
  return(list("mut_vs_wt_in_F" = DEG1,
              "mut_vs_wt_in_M" = DEG2))
  
}
pathway_enrich = function(res){
  geneList = res$t
  names(geneList) = rownames(res)
  geneList = sort(geneList,decreasing = T)
  
  ## gse kegg
  gene_id = bitr(names(geneList), 
                 fromType="SYMBOL", 
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db) 
  names(geneList) = gene_id$ENTREZID
  set.seed(123)
  kk = gseKEGG(geneList     = geneList,
               keyType = "kegg",
               organism     = 'hsa',
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.25,
               verbose      = FALSE,
               seed = 123)
  kk = setReadable(kk,OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  sum(kk@result$p.adjust<0.05)
  sum(kk@result$p.adjust<0.25)
  return(kk)
}

## 1. limma
tcga_node=readRDS("result/tcga_node_score.rds")
TG_all = tcga_node$TG_score_all

tcga_pheno = read.csv("data/TCGA-COADREAD.csv")
tcga_pheno = tcga_pheno[match(colnames(TG_all),tcga_pheno$PATIENT_ID),]
identical(tcga_pheno$PATIENT_ID,colnames(TG_all))

DTGs = tcga_limma(TG_all,tcga_pheno)
saveRDS(DTGs,"result/DTGs_TCGA.rds")

### 2. pathway enrichment
path1 = pathway_enrich(DTGs$mut_vs_wt_in_F)
path2 = pathway_enrich(DTGs$mut_vs_wt_in_M)

write.csv(path1@result,"result/result_final/TCGA_DiffTarg_mut_vs_wt_in_F_kegg_gsea.csv",row.names = F)
write.csv(path2@result,"result/result_final/TCGA_DiffTarg_mut_vs_wt_in_M_kegg_gsea.csv",row.names = F)


## dotplot
gg1 = dotplot(path1,showCategory=11,
              color = "p.adjust",
              x="NES",
              size = NULL,
              title="",
              font.size = 12,
              label_format = 40)
ggsave("plot/plot_final/TCGA_DiffTarg_mut_vs_wt_in_F_kegg_gsea.png",
       gg1,width = 7,height = 3.5,dpi = 600)

## select some representive pathways
## highlighed in TCGA_DiffTarg_mut_vs_wt_in_M_kegg_gsea.xlsx
ID = c("hsa00190",
       "hsa01200",
       "hsa00020",
       "hsa01230",
       "hsa00480",
       "hsa00270",
       "hsa00330",
       "hsa00030",
       "hsa04215",
       "hsa00220",
       "hsa04216"
       )
path2@result = path2@result[path2@result$ID%in%ID,]
gg2 = dotplot(path2,showCategory=11,
             color = "p.adjust",
             x="NES",
             size = NULL,
             title="",
             font.size = 12,
             label_format = 40)
ggsave("plot/plot_final/TCGA_DiffTarg_mut_vs_wt_in_M_kegg_gsea.png",
       gg2,width = 7,height = 3.5,dpi = 600)




