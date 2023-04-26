#****************************************************************************
#*all figures in manuscript by Fig. num
#****************************************************************************

library(TCGAbiolinks)
library(SummarizedExperiment)
#library(sigclust2)
library(dendextend)
library(pheatmap)
library(viridis)
library(survival)
library(survminer)
library(genefu)
library(tidyverse)
library(maftools)
library(ConsensusClusterPlus)
library(AIMS)


# Figure 1 ----------------------------------------------------------------

rm(list = ls())
# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

load("/mnt/Miscrosoft/TCGA_project/Data/R_data/TCGA-BRCA/Expr_and_pheno.RData")
expr_dat <- X1_genename
expr_dat <- log2(expr_dat+1)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]

#find tumor samples
# tumor_sample_ID_dat <- colData(expr) %>% dplyr::as_tibble() %>% 
#   dplyr::filter(sample_type=="Primary Tumor") %>% dplyr::select(barcode,sample_type)
tumor_sample_ID_dat <- A %>% 
  dplyr::filter(sample_type=="Primary Tumor") %>% dplyr::select(barcode,sample_type)


#tumor sample expr
tumor_expr <- expr_dat[,tumor_sample_ID_dat$barcode]


#intrinsic genes
intrinsic_genes <- read.table("2009_JCO_intrinsic_genes_S_table_5.txt",header = FALSE,sep = "\t",
                              stringsAsFactors = FALSE)
#tumor intrinsic gene exprs
#****
#*use the direct intersect between intrinsic genes and genes in dataset,
#*which need to be done more rigorously
#****
tumor_intrinsic_expr <- tumor_expr[intersect(intrinsic_genes[,1],rownames(tumor_expr)),]

# expr_for_clust <- tumor_intrinsic_expr_centering
RNA_seq_expr_for_clust <- tumor_intrinsic_expr


# array data
load("20221216_BRCA_microarray_expression.RData")

expr_dat <- BRCA_array_expr %>% apply(2,as.numeric)
rownames(expr_dat) <- rownames(BRCA_array_expr)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]

#find tumor samples
# tumor_sample_ID_dat <- colData(expr) %>% dplyr::as_tibble() %>% 
#   dplyr::filter(sample_type=="Primary Tumor") %>% dplyr::select(barcode,sample_type)
tumor_sample_ID_dat <- BRCA_array_clin %>% as.data.frame() %>% 
  dplyr::filter(sample_type=="Primary Tumor") %>% dplyr::select(barcode,sample_type)


#tumor sample expr
tumor_expr <- expr_dat[,tumor_sample_ID_dat$barcode]


#intrinsic genes
intrinsic_genes <- read.table("2009_JCO_intrinsic_genes_S_table_5.txt",header = FALSE,sep = "\t",
                              stringsAsFactors = FALSE)
#tumor intrinsic gene exprs
#****
#*use the direct intersect between intrinsic genes and genes in dataset,
#*which need to be done more rigorously
#****
tumor_intrinsic_expr <- tumor_expr[intersect(intrinsic_genes[,1],rownames(tumor_expr)),]

# expr_for_clust <- tumor_intrinsic_expr_centering
array_expr_for_clust <- tumor_intrinsic_expr


#for anno data
wd <- getwd()
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/Basal_like_Breast_Cancer/Results/")
TCGA_symbol_and_entrez <- read.table("TCGA_gene_symbol_and_entrez_ID.txt",
                                     header = TRUE,sep = "\t",stringsAsFactors = FALSE)
setwd(wd)

#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

data("pam50.scale")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50.scale,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


#array data 
tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]


data("pam50.scale")
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50.scale,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_array_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)

clustering_dat_for_plot <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype,
                                 PAM50_subtyping_res_RNA_seq_AIMS$subtype[,1],
                                 PAM50_subtyping_res_array_AIMS$subtype[,1]) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2,RNA_seq_AIMS=3,Array_AIMS=4)


clustering_dat_for_plot %>% select(1,2,5) %>% arrange(Array_PAM50,RNA_seq_PAM50) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank()) +
  labs(x='',y='',fill='Subtypes') 


clustering_dat_for_plot %>% select(3,4,5) %>% arrange(Array_AIMS,RNA_seq_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank()) +
  labs(x='',y='',fill='Subtypes') 



clustering_dat_for_plot %>% arrange(RNA_seq_PAM50,Array_PAM50,RNA_seq_AIMS,Array_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_PAM50","Array_PAM50","RNA_seq_AIMS","Array_AIMS")[4:1])) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank()) +
  labs(x='',y='',fill='Subtypes') 



caret::confusionMatrix(data = factor(clustering_dat_for_plot$RNA_seq_PAM50),
                       reference=factor(clustering_dat_for_plot$Array_PAM50))


caret::confusionMatrix(data = factor(clustering_dat_for_plot$RNA_seq_AIMS),
                       reference=factor(clustering_dat_for_plot$Array_AIMS))

caret::confusionMatrix(data = factor(clustering_dat_for_plot$RNA_seq_PAM50),
                       reference=factor(clustering_dat_for_plot$Array_AIMS))

caret::confusionMatrix(data = factor(clustering_dat_for_plot$RNA_seq_AIMS),
                       reference=factor(clustering_dat_for_plot$Array_PAM50))








# Relative expression model to microarray and RNA-seq ---------------------
#to see the model's consistency
##Using new subtype to subtype TCGA samples

load("../Results/TCGA_train_and_test_intrinsic_gene_AIMS_training_res.RData")


setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")

setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

array_tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

array_pred <- predict.one.vs.all.tsp(D = array_tumor_expr,
                                     GeneName = rownames(array_tumor_expr),
                                     one.vs.all.tsp = train_and_test_res_intrinsic$training_res$final_model)


RNA_seq_tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),
                                                        colnames(array_expr_for_clust))]
RNA_seq_pred <- predict.one.vs.all.tsp(D = RNA_seq_tumor_expr,GeneName = rownames(RNA_seq_tumor_expr),
                                       one.vs.all.tsp = train_and_test_res_intrinsic$training_res$final_model)


clustering_dat_for_plot <- cbind(array_pred$cl,RNA_seq_pred$cl) %>% 
  as.data.frame() %>%  
  rename(Array_New_subtyper=1,RNA_seq_New_subtyper=2) %>% mutate(barcode=rownames(.))

clustering_dat_for_plot %>% arrange(RNA_seq_New_subtyper,Array_New_subtyper) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_New_subtyper","Array_New_subtyper"))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank()) +
  labs(x='',y='',fill='Subtypes') 







