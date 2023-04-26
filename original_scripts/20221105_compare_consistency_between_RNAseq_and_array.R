##classification consistency between RNA-seq and microarray
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sigclust2)
library(tidyverse)
library(dendextend)
library(pheatmap)
library(pamr)
rm(list = ls())

setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")

#using PAM50 
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

load("../Results/TCGA_train_and_test_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))


load("../Results/TCGA_train_and_test_intrinsic_gene_AIMS_training_res.RData")
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))


load("../Results/TCGA_train_and_test_ClaNC_and_PAM_training_res.RData")
#PAM predicting MCC
mean(sapply(PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))
#PAM + Cor predicting MCC
mean(sapply(PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))
#ClaNC predicing
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))
#ClaNC + PAM predicting
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))
#ClaNC + PAM + Cor predicting
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))


# TCGA intrinsic gene ClaNC and PAM ---------------------------------------
load("../Results/TCGA_train_and_test_intrinsic_gene_ClaNC_and_PAM_training_res.RData")
#PAM predicting MCC
mean(sapply(PAM_and_PAM_plus_Cor_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))
#PAM + Cor predicting MCC
mean(sapply(PAM_and_PAM_plus_Cor_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))
#ClaNC predicing
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))
#ClaNC + PAM predicting
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))
#ClaNC + PAM + Cor predicting
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))





# TCGA all gene expr rank MCC ---------------------------------------------

load("../Results/TCGA_train_and_test_by_expr_rank_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))




# TCGA intrinsic gene expr rank MCC ---------------------------------------

load("../Results/TCGA_train_and_test_intrinsic_gene_by_expr_rank_AIMS_training_res.RData")
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))



# TCGA RNA_seq and array training -----------------------------------------

load("../Results/TCGA_RNA_seq_and_array_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))


# TCGA RNA_seq and array training -----------------------------------------

load("../Results/TCGA_RNA_seq_and_array_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))



# TCGA RNA_seq and array intrinsic gene training --------------------------

load("../Results/TCGA_RNA_seq_and_array_intrinsic_AIMS_training_res.RData")
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))








