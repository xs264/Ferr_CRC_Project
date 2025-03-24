library(decoupleR)
library(tidyverse)
library(limma)
library(TIGER)
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
expr = read.csv("data/expression/GSE39582_log2_expr.csv")
expr = column_to_rownames(expr,"geo_accession")
expr = expr[,-c(1:10)]
expr = as.matrix(expr)
expr = t(expr)

# VIPER
net = decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)
TFA_viper = decoupleR::run_viper(mat = expr,
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 1)
TFA_viper = el2adj(TFA_viper[,c(2,3,4)])

# limma
geo_pheno = read.csv("data/GSE39582.csv")
geo_pheno = geo_pheno[match(colnames(TFA_viper),geo_pheno$geo_accession),]
identical(geo_pheno$geo_accession,colnames(TFA_viper))

DTFs = geo_limma(TFA_viper,geo_pheno)

write.csv(DTFs$mut_vs_wt_in_M,"result/result_final/GEO_DiffTFA_mut_vs_wt_in_M.csv")
write.csv(DTFs$mut_vs_wt_in_F,"result/result_final/GEO_DiffTFA_mut_vs_wt_in_F.csv")

## plot
p1 = barplot(DTFs$mut_vs_wt_in_F)
p2 = barplot(DTFs$mut_vs_wt_in_M)

ggsave("plot/plot_final/GEO_DiffTFA_mut_vs_wt_in_F.png",p1,width = 4,height = 4,dpi = 600)
ggsave("plot/plot_final/GEO_DiffTFA_mut_vs_wt_in_M.png",p2,width = 4,height = 4,dpi = 600)









