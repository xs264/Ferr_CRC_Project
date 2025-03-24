library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

geo_limma = function(data,geo_pheno){
  
  mut = geo_pheno$KRAS_mut
  mut[mut=="0: WT"] = "0"
  mut[mut=="1: Mutant"] = "1"
  mut = factor(mut)
  
  sex = factor(geo_pheno$sex)
  age = as.numeric(geo_pheno$age)
  age[which(is.na(age))] = mean(age,na.rm=TRUE)
  stage = factor(geo_pheno$stage)
  chemo = factor(geo_pheno$chemotherapy)
  chemo[is.na(chemo)] = "N"
  
  ### 1. differential targeting
  design = model.matrix(~sex + mut + sex*mut + age + stage + chemo)
  colnames(design) = make.names(colnames(design))
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

### 1. limma
geo_node=readRDS("result/geo_node_score.rds")
TG_all = geo_node$TG_score_all

geo_pheno = read.csv("data/GSE39582.csv")
geo_pheno = geo_pheno[match(colnames(TG_all),geo_pheno$geo_accession),]
identical(geo_pheno$geo_accession,colnames(TG_all))

DTGs = geo_limma(TG_all,geo_pheno)
saveRDS(DTGs,"result/DTGs_GEO.rds")

### 2. pathway enrichment
path1 = pathway_enrich(DTGs$mut_vs_wt_in_F)
path2 = pathway_enrich(DTGs$mut_vs_wt_in_M)

write.csv(path1@result,"result/result_final/GEO_DiffTarg_mut_vs_wt_in_F_kegg_gsea.csv",row.names = F)
write.csv(path2@result,"result/result_final/GEO_DiffTarg_mut_vs_wt_in_M_kegg_gsea.csv",row.names = F)

### 3. dotplot
gg1 = dotplot(path1,showCategory=11,
              color = "p.adjust",
              x="NES",
              size = NULL,
              title="",
              font.size = 12,
              label_format = 40)
ggsave("plot/plot_final/GEO_DiffTarg_mut_vs_wt_in_F_kegg_gsea.png",
       gg1,width = 7,height = 3.5,dpi = 600)

gg2 = dotplot(path2,showCategory=11,
              color = "p.adjust",
              x="NES",
              size = NULL,
              title="",
              font.size = 12,
              label_format = 40)
ggsave("plot/plot_final/GEO_DiffTarg_mut_vs_wt_in_M_kegg_gsea.png",
       gg2,width = 7,height = 3.5,dpi = 600)








