#****************************************************************************
#*Supplementary Figure 3 in manuscript
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
#library(WGCNA)
library(dendextend)


# Figure 1 ----------------------------------------------------------------

rm(list = ls())
# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("./Data/")

load("./R_data/BRCA/Expr_and_pheno.RData")
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













# For clustering intrinsic gene clustering --------------------------------

RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[,
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,scale))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,scale))

combined_dat <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_RNAseq"),
                            paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_array"))

clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D") %>% as.dendrogram()





#***************************
#*Supplementary Figure 2****
#***************************


pdf("../Results/Figure/Supplementary_Figure_3a.pdf",width = 20,height = 10)
clustering_res %>%  
  dendextend::set("labels_colors",ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8")) %>%
  #set("branches_col",ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8")) %>% 
  assign_values_to_leaves_edgePar(value = ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8"),
                                  edgePar = "col") %>%
  set("labels_cex",0.1) %>% 
  plot(main="Z-Scaled data",cex.main=2)
legend("topright",legend = c("Microarray","RNA-seq"),col=c("#E41A1C","#377EB8"),lwd=1,bty = "n",cex = 3)
dev.off()





pdf("../Results/Figure/Supplementary_Figure_3b.pdf",width = 20,height = 10)
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,function(x)x-median(x)))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,function(x)x-median(x)))

combined_dat1 <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat1) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                    "_RNAseq"),
                             paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                    "_array"))

clustering_res <- hclust(d = dist(t(combined_dat1)),method = "ward.D") %>% as.dendrogram()
clustering_res %>%  
  dendextend::set("labels_colors",ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8")) %>%
  #set("branches_col",ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8")) %>% 
  assign_values_to_leaves_edgePar(value = ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8"),
                                  edgePar = "col") %>%
  set("labels_cex",0.1) %>% 
  plot(main="Median-Centered data",cex.main=2)
legend("topright",legend = c("Microarray","RNA-seq"),col=c("#E41A1C","#377EB8"),lwd=1,bty = "n",cex = 3)

dev.off()


pdf("../Results/Figure/Supplementary_Figure_3c.pdf",width = 20,height = 10)
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,function(x)x-mean(x)))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,function(x)x-mean(x)))

combined_dat1 <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat1) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                    "_RNAseq"),
                             paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                    "_array"))
clustering_res <- hclust(d = dist(t(combined_dat1)),method = "ward.D") %>% as.dendrogram()
clustering_res %>%  
  dendextend::set("labels_colors",ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8")) %>%
  #set("branches_col",ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8")) %>% 
  assign_values_to_leaves_edgePar(value = ifelse(grepl("_array",labels(clustering_res)),"#E41A1C","#377EB8"),
                                  edgePar = "col") %>%
  set("labels_cex",0.1) %>% 
  plot(main="Mean-Centered data",cex.main=2)
legend("topright",legend = c("Microarray","RNA-seq"),col=c("#E41A1C","#377EB8"),lwd=1,bty = "n",cex = 3)
dev.off()



















