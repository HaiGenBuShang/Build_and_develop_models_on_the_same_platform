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
library(WGCNA)


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



# For clustering PAM50 genes ----------------------------------------------

data("pam50")
PAM50_genes <- rownames(pam50$centroids)

RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[intersect(PAM50_genes,rownames(RNA_seq_expr_for_clust)),
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[intersect(PAM50_genes,rownames(array_expr_for_clust)),
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]


intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))


RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,scale))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,scale))

# RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,function(x)x-mean(x)))
# array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,function(x)x-mean(x)))


combined_dat <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_RNAseq"),
                            paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_array"))

annotation_col <- data.frame(data_type=rep(c("RNA-seq","Array"),
                                           times=c(ncol(RNA_seq_expr_PAM50),ncol(array_expr_PAM50))),
                             samples=intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)))
rownames(annotation_col) <- colnames(combined_dat)

ann_colors <- list(data_type= c(`RNA-seq` ="#A04EF6", Array = "#E7D044"),
                   samples=setNames(rainbow(ncol(RNA_seq_expr_PAM50)),
                                    intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))))
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),"TCGA-",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "complete",
         annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))


annotation_col <- data.frame(data_type=rep(c("RNA-seq","Array"),
                                           times=c(ncol(RNA_seq_expr_PAM50),ncol(array_expr_PAM50))))
rownames(annotation_col) <- colnames(combined_dat)



pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),".*",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "complete",
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))


clustering_res <- hclust(d = dist(t(combined_dat)),method = "complete")

samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1)


samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n!=2)








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

annotation_col <- data.frame(data_type=rep(c("RNA-seq","Array"),
                                           times=c(ncol(RNA_seq_expr_PAM50),ncol(array_expr_PAM50))),
                             samples=intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)))
rownames(annotation_col) <- colnames(combined_dat)

ann_colors <- list(data_type= c(`RNA-seq` ="#A04EF6", Array = "#E7D044"),
                   samples=setNames(rainbow(ncol(RNA_seq_expr_PAM50)),
                                    intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))))
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),"TCGA-",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "complete",
         annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))


annotation_col <- data.frame(data_type=rep(c("RNA-seq","Array"),
                                           times=c(ncol(RNA_seq_expr_PAM50),ncol(array_expr_PAM50))))
rownames(annotation_col) <- colnames(combined_dat)
ann_colors <- list(data_type= c(`RNA-seq` ="#A04EF6", Array = "#E7D044"))


pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),".*",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "complete",annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))


clustering_res <- hclust(d = dist(t(combined_dat)),method = "complete")

samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1)


samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n!=2)





# For original expr intrinsic gene clustering -----------------------------

RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[,
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))
# RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,scale))
# array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,scale))

RNA_seq_expr_PAM50 <- RNA_seq_expr_PAM50[intersect_genes,]
array_expr_PAM50 <- array_expr_PAM50[intersect_genes,]

combined_dat <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_RNAseq"),
                            paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_array"))

annotation_col <- data.frame(data_type=rep(c("RNA-seq","Array"),
                                           times=c(ncol(RNA_seq_expr_PAM50),ncol(array_expr_PAM50))),
                             samples=intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)))
rownames(annotation_col) <- colnames(combined_dat)

ann_colors <- list(data_type= c(`RNA-seq` ="#A04EF6", Array = "#E7D044"),
                   samples=setNames(rainbow(ncol(RNA_seq_expr_PAM50)),
                                    intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))))
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),"TCGA-",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "complete",
         annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))


annotation_col <- data.frame(data_type=rep(c("RNA-seq","Array"),
                                           times=c(ncol(RNA_seq_expr_PAM50),ncol(array_expr_PAM50))))
rownames(annotation_col) <- colnames(combined_dat)
ann_colors <- list(data_type= c(`RNA-seq` ="#A04EF6", Array = "#E7D044"))


pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),".*",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "complete",annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))


clustering_res <- hclust(d = dist(t(combined_dat)),method = "complete")

samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1)


samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n!=2)








#gene scale for AIMS 
RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[,
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,function(x)x-mean(x)))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,function(x)x-mean(x)))
combined_dat <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_RNAseq"),
                            paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_array"))
#for anno data
wd <- getwd()
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/Basal_like_Breast_Cancer/Results/")
TCGA_symbol_and_entrez <- read.table("TCGA_gene_symbol_and_entrez_ID.txt",
                                     header = TRUE,sep = "\t",stringsAsFactors = FALSE)
setwd(wd)

tumor_expr <- 2^combined_dat
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


PAM50_subtyping_res_RNA_seq_AIMS$subtype %>% as.data.frame %>% mutate(barcode=rownames(.)) %>% 
  mutate(samples=str_replace_all(barcode,"_.*","")) %>% group_by(samples,"20") %>% count() %>% 
  filter(n==2)





# *************************************************************************
# After Scale, the prediction are similiar --------------------------------
# *************************************************************************




#Original value subtyping using PAM50
RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[,
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

#RNA_seq data
tumor_expr <- RNA_seq_expr_PAM50
data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)
#array data
tumor_expr <- array_expr_PAM50
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)


cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype) %>% as.data.frame() %>% 
  mutate(barcode=rownames(.))  %>% mutate(same_assign=V1==V2) %>% count(same_assign)




















