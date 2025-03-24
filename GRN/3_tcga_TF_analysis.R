library(TIGER)
library(dorothea)
library(decoupleR)
library(tidyverse)
library(limma)
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
barplot = function(contrast_acts){
  
  contrast_acts$"statistics" = "viper"
  contrast_acts$"source" = rownames(contrast_acts)
  contrast_acts$"condition" = "t"
  contrast_acts$"score" = contrast_acts$t
  contrast_acts$"p_value" = contrast_acts$P.Value
  
  f_contrast_acts <- contrast_acts %>%
    dplyr::mutate(rnk = NA)
  
  msk <- f_contrast_acts$score > 0
  
  f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
  f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
  
  tfs <- f_contrast_acts %>%
    dplyr::arrange(rnk) %>%
    head(10) %>%
    dplyr::pull(source)
  
  f_contrast_acts <- f_contrast_acts %>%
    filter(source %in% tfs)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  
  p <- ggplot2::ggplot(data = f_contrast_acts, 
                       mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                              y = score)) + 
    ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                      color = "black",
                      stat = "identity") +
    ggplot2::scale_fill_gradient2(low = colors[1], 
                                  mid = "whitesmoke", 
                                  high = colors[2], 
                                  midpoint = 0) + 
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                   axis.text.x = ggplot2::element_text(angle = 45, 
                                                       hjust = 1, 
                                                       size = 10, 
                                                       face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 10, 
                                                       face = "bold"),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank()) +
    ggplot2::xlab("TFs")
  
  
  return(p)
}

# data
expr = read.csv("data/expression/TCGA-COADREAD_log2_expr.csv")
expr = column_to_rownames(expr,"PATIENT_ID")
expr = expr[,-c(1:8)]
expr = as.matrix(expr)
expr = t(expr)


# VIPER
net = decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)
TFA_viper = decoupleR::run_viper(#mat = DTGs$FvsM_in_mut1[,"t",drop=F],
                                  mat = expr,
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 1)
TFA_viper = el2adj(TFA_viper[,c(2,3,4)])


# limma
tcga_pheno = read.csv("data/TCGA-COADREAD.csv")
tcga_pheno = tcga_pheno[match(colnames(TFA_viper),tcga_pheno$PATIENT_ID),]
identical(tcga_pheno$PATIENT_ID,colnames(TFA_viper))

DTFs = tcga_limma(TFA_viper,tcga_pheno)
saveRDS(DTFs,"result/DTF_activity_TCGA.rds")
write.csv(DTFs$mut_vs_wt_in_F,"result/result_final/TCGA_DiffTFA_mut_vs_wt_in_F.csv")
write.csv(DTFs$mut_vs_wt_in_M,"result/result_final/TCGA_DiffTFA_mut_vs_wt_in_M.csv")

## plot
p1 = barplot(DTFs$mut_vs_wt_in_F)
p2 = barplot(DTFs$mut_vs_wt_in_M)

ggsave("plot/plot_final/TCGA_DiffTFA_mut_vs_wt_in_F.png",p1,width = 4,height = 4,dpi = 600)
ggsave("plot/plot_final/TCGA_DiffTFA_mut_vs_wt_in_M.png",p2,width = 4,height = 4,dpi = 600)




