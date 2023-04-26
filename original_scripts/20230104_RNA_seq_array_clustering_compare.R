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
library(AIMS)

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

RNA_seq_expr_rank_for_clust <- apply(tumor_intrinsic_expr,2,function(x){
  rank(-1*x)
})



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
array_expr_rank_for_clust <- apply(tumor_intrinsic_expr,2,function(x){
  rank(-1*x)
})




# RNA_seq data clustering using the sample samples in array ---------------
expr_for_clust <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

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
#rm(list = grep("RNA_seq_clustering_res",ls(),value = TRUE,invert = TRUE))


# array data clustering using the sample samples in RNA_seq ---------------
expr_for_clust <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

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

clustering_res <- array_clustering_res %>% mutate(barcode=rownames(.)) %>% 
  inner_join(RNA_seq_clustering_res %>% mutate(barcode=rownames(.)),by="barcode",suffix=c(".array",".RNA_seq"))

# compare_two_clusterings <- function(clustering_dat,clustering_col_1,clustering_col_2,
#                                     if_return_info=FALSE){
#   clustering_1 <- clustering_dat %>% group_by_at(clustering_col_1) %>% 
#     summarise(clustering_res=paste(sort(barcode),collapse = ", ")) %>% 
#     mutate(samples=map(clustering_res,~(str_split(.,pattern = ", "))[[1]]))
#   clustering_2 <- clustering_dat %>% group_by_at(clustering_col_2) %>% 
#     summarise(clustering_res=paste(sort(barcode),collapse = ", ")) %>% 
#     mutate(samples=map(clustering_res,~(str_split(.,pattern = ", "))[[1]]))
#   
#   matching_dat <- sapply(clustering_1$samples,function(x){
#     sapply(clustering_2$samples,function(y){
#       c(length(intersect(x,y))/length(x),length(intersect(x,y))/length(y))
#     })
#   },simplify = FALSE)
#   
#   match_results <- matching_dat %>% lapply(function(x){
#     c(which.max(x[1,]),which.max(x[2,]))
#   })
#   
#   invisible(sapply(1:length(match_results),function(x){
#     cat(paste0("The ",x," cluster in clustering 1 corresponding to the ",
#                unique(match_results[[x]])," clustering in clustering 2\n"))
#   }))
#   if(if_return_info){
#     matching_dat
#   }else{
#     invisible(matching_dat)
#   }
# }
# 
# compare_two_clusterings(clustering_dat = clustering_res,
#                         clustering_col_1 = "cluster_pearson.array",clustering_col_2 = "clutser_euclidean.RNA_seq")
# clustering_dat_for_plot <- clustering_res %>%
#   mutate(clutser_euclidean.RNA_seq = as.integer(factor(clutser_euclidean.RNA_seq,levels = c(1,3,2,5,4))))

(clustering_res[,c("cluster_pearson.array","clutser_euclidean.RNA_seq")] %>% 
    group_by(cluster_pearson.array) %>% nest() %>% 
    mutate(x=map(data,~(table(.)))))$x #different clusters mutually contained
clustering_dat_for_plot <- clustering_res %>% 
  mutate(clutser_euclidean.RNA_seq = as.integer(factor(clutser_euclidean.RNA_seq,levels = c(5,3,4,2,1))))

# clustering_dat_for_plot <- clustering_res %>% 
#   mutate(clutser_euclidean.RNA_seq = as.integer(factor(clutser_euclidean.RNA_seq)))



# clustering_dat_for_plot %>% select(3:5) %>% 
#   pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
#   arrange(Clustering,value) %>% 
#   mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
#   ggplot() +
#   geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
#   #scale_fill_viridis(discrete=FALSE) +
#   scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Set3"),"red","green","blue")) +
#   theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
#   labs(x='',y='',fill='Coef in ARD')


clustering_dat_for_plot %>% select(3:5) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(Clustering,value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
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






# AIMS predicting results for RNA_seq and array data ----------------------

#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

data("pam50.robust")
#note that different normalization method would affect the results
PAM50_subtyping_res_RNA_seq <- molecular.subtyping(sbt.model = "pam50",data = t(tumor_expr),
                                                   annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                      stringsAsFactors = FALSE),
                                                   do.mapping = FALSE)
data("pam50.scale")
PAM50_subtyping_res_RNA_seq <- intrinsic.cluster.predict(sbt.model = pam50.scale,data = t(tumor_expr),
                                                         annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                            stringsAsFactors = FALSE),
                                                         do.mapping = FALSE)

#array data
tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
tumor_expr <- 2^tumor_expr

#for anno data
wd <- getwd()
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/Basal_like_Breast_Cancer/Results/")
TCGA_symbol_and_entrez <- read.table("TCGA_gene_symbol_and_entrez_ID.txt",
                                     header = TRUE,sep = "\t",stringsAsFactors = FALSE)
setwd(wd)
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]

anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]

PAM50_subtyping_res_array <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                 annot = anno,
                                                 do.mapping = FALSE)
all(rownames(PAM50_subtyping_res_RNA_seq$subtype %>% as.data.frame())==rownames(PAM50_subtyping_res_array$subtype))

cbind(PAM50_subtyping_res_array$subtype, PAM50_subtyping_res_RNA_seq$subtype %>% as.data.frame())

table(PAM50_subtyping_res_array$subtype[,1],PAM50_subtyping_res_RNA_seq$subtype)
caret::confusionMatrix(data=factor(PAM50_subtyping_res_array$subtype[,1]),
                       reference=factor(PAM50_subtyping_res_RNA_seq$subtype))









# RNA_seq clustering using expr rank --------------------------------------

expr_for_clust <- RNA_seq_expr_rank_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),
                                                         colnames(array_expr_for_clust))]

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
#rm(list = grep("RNA_seq_clustering_res",ls(),value = TRUE,invert = TRUE))


# array clustering using expr rank ----------------------------------------

expr_for_clust <- array_expr_rank_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),
                                                       colnames(array_expr_for_clust))]

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


clustering_res <- array_clustering_res %>% mutate(barcode=rownames(.)) %>% 
  inner_join(RNA_seq_clustering_res %>% mutate(barcode=rownames(.)),by="barcode",suffix=c(".array",".RNA_seq"))


(clustering_res[,c("cluster_pearson.array","clutser_euclidean.RNA_seq")] %>% 
    group_by(cluster_pearson.array) %>% nest() %>% 
    mutate(x=map(data,~(table(.)))))$x #different clusters mutually contained


clustering_dat_for_plot <- clustering_res %>% 
  mutate(clutser_euclidean.RNA_seq = as.integer(factor(clutser_euclidean.RNA_seq,levels = c(1,2,3,5,4))))

clustering_dat_for_plot %>% select(3:5) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(Clustering,value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='',fill='Coef in ARD')

caret::confusionMatrix(data=factor(clustering_dat_for_plot$cluster_pearson.array),
                       reference=factor(clustering_dat_for_plot$clutser_euclidean.RNA_seq))



clustering_dat_for_plot <- clustering_res %>% 
  mutate(clutser_euclidean.RNA_seq = as.integer(factor(clutser_euclidean.RNA_seq,levels = c(2,3,4,5,1))))

clustering_dat_for_plot %>% select(1,4,5) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>%
  arrange(Clustering,value) %>% 
  mutate(barcode=factor(barcode,factor(unique(barcode)))) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='',fill='Coef in ARD')

caret::confusionMatrix(data=factor(clustering_dat_for_plot$clutser_euclidean.array),
                       reference=factor(clustering_dat_for_plot$clutser_euclidean.RNA_seq))
