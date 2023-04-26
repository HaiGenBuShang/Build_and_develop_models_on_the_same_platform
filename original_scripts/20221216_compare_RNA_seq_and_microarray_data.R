#****************************************************************************
#*to compare clustering results between RNA-seq and microarray
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

set.seed(12345678)
#down load TCGA BRCA expressions
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

# RNA-seq data clustering -------------------------------------------------

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
RNA_seq_tumor_expr <- tumor_expr

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
expr_for_clust <- tumor_intrinsic_expr

#do clust
#use pearson and average
cors <- cor(expr_for_clust,use = "pairwise.complete.obs", method = "pearson")
clust_res <- hclust(d = as.dist(1-cors),method = "average")

#cut to 5 classes
samples_of_classes <- cutree(clust_res,k = 5)
#plot(color_branches(clust_res,k=5),leaflab="none")

gene_clust <- hclust(d = dist(expr_for_clust,method = "euclidean"),method = "ward.D")

#do clust
clust_res_euclidean <- hclust(d = dist(t(expr_for_clust),method = "euclidean"),method = "ward.D")
samples_of_classes_euclidean <- cutree(clust_res_euclidean,k=5)

#do clust by median centering
clust_dat_centering <- t(apply(expr_for_clust,1,function(x){
  x-median(x)
}))

cors_centering <- cor(clust_dat_centering,use = "pairwise.complete.obs", method = "pearson")
clust_res_centering <- hclust(d = as.dist(1-cors_centering),method = "average")
samples_of_classes_centering <- cutree(clust_res_centering,k=5)

RNA_seq_clustering_res <- data.frame(clutser_euclidean=samples_of_classes_euclidean,
                                     cluster_pearson_centering=samples_of_classes,
                                     cluster_pearson=samples_of_classes_centering)
rm(list = grep("RNA_seq_clustering_res",ls(),value = TRUE,invert = TRUE))




# Microarray data clustering ----------------------------------------------

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
array_tumor_expr <- tumor_expr

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
expr_for_clust <- tumor_intrinsic_expr

#do clust
#use pearson and average
cors <- cor(expr_for_clust,use = "pairwise.complete.obs", method = "pearson")
clust_res <- hclust(d = as.dist(1-cors),method = "average")

#cut to 5 classes
samples_of_classes <- cutree(clust_res,k = 5)
#plot(color_branches(clust_res,k=5),leaflab="none")

gene_clust <- hclust(d = dist(expr_for_clust,method = "euclidean"),method = "ward.D")

#do clust
clust_res_euclidean <- hclust(d = dist(t(expr_for_clust),method = "euclidean"),method = "ward.D")
samples_of_classes_euclidean <- cutree(clust_res_euclidean,k=5)

#do clust by median centering
clust_dat_centering <- t(apply(expr_for_clust,1,function(x){
  x-median(x)
}))

cors_centering <- cor(clust_dat_centering,use = "pairwise.complete.obs", method = "pearson")
clust_res_centering <- hclust(d = as.dist(1-cors_centering),method = "average")
samples_of_classes_centering <- cutree(clust_res_centering,k=5)

array_clustering_res <- data.frame(clutser_euclidean=samples_of_classes_euclidean,
                                   cluster_pearson_centering=samples_of_classes,
                                   cluster_pearson=samples_of_classes_centering,
                                   row.names = names(samples_of_classes))

#rm(list = grep("RNA_seq_clustering_res|array_clustering_res",ls(),value = TRUE,invert = TRUE))

clustering_res <- array_clustering_res %>% mutate(barcode=rownames(.)) %>% 
  inner_join(RNA_seq_clustering_res %>% mutate(barcode=rownames(.)),by="barcode",suffix=c(".array",".RNA_seq"))

compare_two_clusterings <- function(clustering_dat,clustering_col_1,clustering_col_2,
                                    if_return_info=FALSE){
  clustering_1 <- clustering_dat %>% group_by_at(clustering_col_1) %>% 
    summarise(clustering_res=paste(sort(barcode),collapse = ", ")) %>% 
    mutate(samples=map(clustering_res,~(str_split(.,pattern = ", "))[[1]]))
  clustering_2 <- clustering_dat %>% group_by_at(clustering_col_2) %>% 
    summarise(clustering_res=paste(sort(barcode),collapse = ", ")) %>% 
    mutate(samples=map(clustering_res,~(str_split(.,pattern = ", "))[[1]]))
  
  matching_dat <- sapply(clustering_1$samples,function(x){
    sapply(clustering_2$samples,function(y){
      c(length(intersect(x,y))/length(x),length(intersect(x,y))/length(y))
    })
  },simplify = FALSE)
  
  match_results <- matching_dat %>% lapply(function(x){
    c(which.max(x[1,]),which.max(x[2,]))
  })
  
  invisible(sapply(1:length(match_results),function(x){
    cat(paste0("The ",x," cluster in clustering 1 corresponding to the ",
               unique(match_results[[x]])," clustering in clustering 2\n"))
  }))
  if(if_return_info){
    matching_dat
  }else{
    invisible(matching_dat)
  }
}

compare_two_clusterings(clustering_dat = clustering_res,
                        clustering_col_1 = "cluster_pearson.array",clustering_col_2 = "clutser_euclidean.RNA_seq")
clustering_dat_for_plot <- clustering_res %>% 
  mutate(clutser_euclidean.RNA_seq = as.integer(factor(clutser_euclidean.RNA_seq,levels = c(1,3,2,5,4))))




clustering_dat_for_plot %>% select(3:5) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(Clustering,value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Set3"),"red","green","blue")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='',fill='Coef in ARD')



clustering_res %>% select(3,4,7) %>% pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(Clustering,value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Set3"),"red","green","blue")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='')






# Relative expression model to microarray and RNA-seq ---------------------
#to see the model's consistency
##Using new subtype to subtype TCGA samples

#load model
#load("../Results/TCGA_train_and_testAIMS_training_res.RData")

load("../Results/TCGA_train_and_test_intrinsic_gene_AIMS_training_res.RData")


setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")

# setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
# array_pred <- predict.one.vs.all.tsp(D = tumor_expr,GeneName = rownames(tumor_expr),
#                                      one.vs.all.tsp = AIMS_training_res$training_res$final_model)

setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
array_pred <- predict.one.vs.all.tsp(D = array_tumor_expr,GeneName = rownames(array_tumor_expr),
                                     one.vs.all.tsp = train_and_test_res_intrinsic$training_res$final_model)
array_pred$cl


# RNA_seq_pred <- predict.one.vs.all.tsp(D = tumor_expr,GeneName = rownames(tumor_expr),
#                                        one.vs.all.tsp = AIMS_training_res$training_res$final_model)

RNA_seq_pred <- predict.one.vs.all.tsp(D = RNA_seq_tumor_expr,GeneName = rownames(RNA_seq_tumor_expr),
                                       one.vs.all.tsp = train_and_test_res_intrinsic$training_res$final_model)

array_pred$cl %>% data.frame() %>% mutate(barcode=rownames(.)) %>% head()
RNA_seq_pred$cl %>% data.frame() %>% mutate(barcode=rownames(.)) %>% head()

relative_expr_predicting <- array_pred$cl %>% data.frame() %>% mutate(barcode=rownames(.)) %>% 
  inner_join(RNA_seq_pred$cl %>% data.frame() %>% mutate(barcode=rownames(.)),by="barcode",
             suffix=c(".array",".RNA_seq"))

### do not need to guess because they were predicted
# compare_two_clusterings(clustering_dat = relative_expr_predicting,
#                         clustering_col_1 = "X46.array",clustering_col_2 = "X46.RNA_seq")

clustering_dat_for_plot <- relative_expr_predicting %>% 
  mutate(X46.array = as.integer(X46.array)) %>% 
  mutate(X46.RNA_seq = as.integer(X46.RNA_seq))

clustering_dat_for_plot %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(Clustering,value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  #scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Set3"),"red","green","blue")) +
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='',fill='Coef in ARD')


clustering_dat_for_plot %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(desc(Clustering),value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  #scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Set3"),"red","green","blue")) +
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='',fill='Coef in ARD')


caret::confusionMatrix(data=factor(clustering_dat_for_plot$X46.array),
                       reference=factor(clustering_dat_for_plot$X46.RNA_seq))





