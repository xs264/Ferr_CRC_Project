
setwd("C:/Snow/CCM draft/CCM draft/VIMP results")
data<-read.csv("vimpCI_male_tcga_all_age.csv",stringsAsFactors = T, header = T)
colnames(data)[1]<-'gene'
data<-as.data.frame.matrix(data)

library(dplyr)
target<-c("FTH1","NFS1","FH","GPX4","ACSL3","FTL","SLC38A1","PPP1R13L", "NR4A1")
genes <- filter(data, gene%in%target)

data1<-read.csv("vimpCI_female_tcga_all_age.csv",stringsAsFactors = T, header = T)
colnames(data1)[1]<-'gene'
data1<-as.data.frame.matrix(data1)
genes1 <- filter(data1, gene%in%target)

gene<-merge(genes,genes1,by="gene")

library(ggplot2)
p1<-genes %>%
  arrange(pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(name=factor(gene, levels=gene)) %>%   # This trick update the factor levels
  ggplot(aes(x=gene, y=-log10(pvalue)))+
  geom_segment( aes(xend=gene, yend=0)) +
  geom_point( size=4, color="#E5ADCB") +
  coord_flip() +
  theme_bw() +
  xlab("Suppressor genes")+ylab("Permutation -log10(P value)")+
  theme(text = element_text(size = 16))+ggtitle("Male")+
  geom_hline(yintercept=1.30103, linetype="dashed", color = "red", size=1.5)

p2<-genes1 %>%
  arrange(pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(name=factor(gene, levels=gene)) %>%   # This trick update the factor levels
  ggplot( aes(x=gene, y=-log10(pvalue))) +
  geom_segment( aes(xend=gene, yend=0)) +
  geom_point( size=4, color="#955677") +
  coord_flip() +
  theme_bw() +
  xlab("Suppressor genes")+ylab("Permutation -log10(P value)")+
  theme(text = element_text(size = 16))+ggtitle("Female")+
  geom_hline(yintercept=1.30103, linetype="dashed", color = "red", size=1.5)

## driver
data2<-read.csv("vimpCI_male_tcga_all_age.csv",stringsAsFactors = T, header = T)
colnames(data2)[1]<-'gene'
data2<-as.data.frame.matrix(data2)

library(dplyr)
target<-c("SLC1A5","ELOVL5","ALOX5","ALOX12","ALOX12B","ALOX15","ACSL4")
genes2 <- filter(data2, gene%in%target)

data3<-read.csv("vimpCI_female_tcga_all_age.csv",stringsAsFactors = T, header = T)
colnames(data3)[1]<-'gene'
data3<-as.data.frame.matrix(data3)
genes3 <- filter(data3, gene%in%target)

gene<-merge(genes2,genes3,by="gene")

library(ggplot2)
p3<-genes2 %>%
  arrange(pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(name=factor(gene, levels=gene)) %>%   # This trick update the factor levels
  ggplot( aes(x=gene, y=-log10(pvalue))) +
  geom_segment( aes(xend=gene, yend=0)) +
  geom_point( size=4, color="#E5ADCB") +
  coord_flip() +
  theme_bw() +
  xlab("Driver genes")+ylab("Permutation -log10(P value)")+
  theme(text = element_text(size = 16))+ggtitle("Male")+
  geom_hline(yintercept=1.30103, linetype="dashed", color = "red", size=1.5)

p4<-genes3 %>%
  arrange(pvalue) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(name=factor(gene, levels=gene)) %>%   # This trick update the factor levels
  ggplot( aes(x=gene, y=-log10(pvalue))) +
  geom_segment( aes(xend=gene, yend=0)) +
  geom_point( size=4, color="#955677") +
  coord_flip() +
  theme_bw() +
  xlab("Driver genes")+ylab("Permutation -log10(P value)")+
  theme(text = element_text(size = 16))+ggtitle("Female")+
  geom_hline(yintercept=1.30103, linetype="dashed", color = "red", size=1.5)

library(ggpubr)
p<-ggarrange(p1,p2,p3,p4,#labels = c("A", "B","C","D"),
             nrow=2, ncol=2,common.legend = TRUE, legend = "bottom")

ggsave(
  filename = ("barplot_TCGA.jpg"),
  plot = p,
  width = 12,
  height = 10,
  dpi=600,
  device = "jpg"
)

