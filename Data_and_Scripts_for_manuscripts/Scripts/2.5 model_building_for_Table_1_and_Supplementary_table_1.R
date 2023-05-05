#****************************************************************************
#*build model for supplementary Table 1
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



message(paste0("Make sure the package \"e107\" \"gplots\" \"ROCR\" and \"multicore\" were installed!\n",
               "And \"Rgtsp\", which located in the AIMS package, were also installed!\n",
               "Besides, make sure there are no duplcated rownames in the expression matrix!!!"))


setwd("../Scripts/")
source("./clanc.R")
source("./all_functions_20220302.R")
source("./trainAIMS_2.R")


#*****************
#*Define functions
#*****************
AIMS_train_and_pred <- function(train_and_test_data_with_lable,PREFIX,k.fold,num.of.rules,...){
  trained_model <- trainAIMS_2(D = train_and_test_data_with_lable$train_set$x,
                               cl = as.character(train_and_test_data_with_lable$train_set$y),
                               EntrezID = rownames(train_and_test_data_with_lable$train_set$x),PREFIX = PREFIX,
                               k.fold = k.fold,num.of.rules = num.of.rules,...)
  predicted_res <- predict.one.vs.all.tsp(D = train_and_test_data_with_lable$test_set$x,
                                          GeneName = rownames(train_and_test_data_with_lable$train_set$x),
                                          one.vs.all.tsp = trained_model$final_model)
  test_set_confusion <- cbind(predicted_res$cl,as.character(train_and_test_data_with_lable$test_set$y))
  colnames(test_set_confusion) <- c("predict","real")
  AIMS_training_res <- list(training_res=trained_model,predicted_res=predicted_res,
                            test_confusion=table(as.data.frame(test_set_confusion,stringsAsFactors=FALSE)),
                            train_data=train_and_test_data_with_lable$train_set$x)
  #save(AIMS_training_res,file = paste0(PREFIX,"AIMS_training_res.RData"))
  AIMS_training_res
}


#rm(list = ls())
# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("../Data/")

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
#RNA_seq_expr_for_clust <- tumor_intrinsic_expr
RNA_seq_expr_for_clust <- tumor_expr


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
#array_expr_for_clust <- tumor_intrinsic_expr
array_expr_for_clust <- tumor_expr


#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
# tumor_expr <- tumor_expr[apply(tumor_expr,1,function(x){
#   !any(is.na(x))
# }),]
# tumor_expr <- tumor_expr[apply(tumor_expr,1,sd)!=0,]
#scale RNA-seq data
tumor_expr <- t(apply(tumor_expr,1,scale))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)


#array data 
tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
# tumor_expr <- tumor_expr[apply(tumor_expr,1,function(x){
#   !any(is.na(x))
# }),]
# tumor_expr <- tumor_expr[apply(tumor_expr,1,sd)!=0,]

#scale array data
tumor_expr <- t(apply(tumor_expr,1,scale))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                             annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                stringsAsFactors = FALSE),
                                                             do.mapping = FALSE)


clustering_dat_for_plot <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2)



RNA_seq_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
intersect_genes <- intersect(rownames(RNA_seq_expr),rownames(array_expr))



tumor_expr <- cbind(RNA_seq_expr[intersect_genes,],array_expr[intersect_genes,])
colnames(tumor_expr) <- c(paste0("RNA_seq_",
                                 intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))),
                          paste0("Array_",
                                 intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))))
tumor_expr <- tumor_expr[apply(tumor_expr,1,function(x){
  !any(is.na(x))
}),]
tumor_expr <- tumor_expr[apply(tumor_expr,1,sd)!=0,]


cluster_sample_subtypes <- c(clustering_dat_for_plot$RNA_seq_PAM50,clustering_dat_for_plot$Array_PAM50)


#train model using AIMS idea
dat_for_model <- tumor_expr
datset_for_model <- list(x=dat_for_model,y=cluster_sample_subtypes)
datset_for_model_intrinsic_g <- list(x=dat_for_model[intersect(intrinsic_genes[,1],rownames(dat_for_model)),],
                                     y=cluster_sample_subtypes)
res_prefix <- "TCGA_RNA_seq_and_array"
res_prefix_intrinsic <- "TCGA_RNA_seq_and_array_intrinsic"




#using all genes
setwd("../Results/")

#****************
#****NOTE!!!!****
#****NOTE!!!!****
#****NOTE!!!!****
#****************

#***********************************************************************************
#*The different scripts must be run in different clean R environment
#*otherwise it will cost lots of time or it would throw errors out
#*You should finishing run one script part, e.g. Part 1, then quit R
#*and then run the commands above and skip Part 1, and run Part 2, and then quit R
#*Until you finish runing all the script parts
#***********************************************************************************




#**************
#*Script Part 1
#**************
set.seed(12345678)
train_and_test_with_label <- produce_train_test_set(expr_with_label = datset_for_model)
train_and_test_res <- AIMS_train_and_pred(train_and_test_data_with_lable = train_and_test_with_label,
                                          PREFIX = res_prefix,
                                          k.fold = 10,num.of.rules = seq(1,50,1))
save(train_and_test_res,file = paste0(res_prefix,"_AIMS_training_res.RData"))
message("Finished all genes AIMS!")

# #**************
# #*Script Part 2
# #**************
# #using intrinsic genes
# train_and_test_with_label_intrinsic <- produce_train_test_set(expr_with_label = datset_for_model_intrinsic_g)
# train_and_test_res_intrinsic<-AIMS_train_and_pred(
#   train_and_test_data_with_lable=train_and_test_with_label_intrinsic,PREFIX = res_prefix_intrinsic,
#   k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res_intrinsic,file = paste0(res_prefix_intrinsic,"_AIMS_training_res.RData"))
# message("Finished intrinsic genes AIMS!")





