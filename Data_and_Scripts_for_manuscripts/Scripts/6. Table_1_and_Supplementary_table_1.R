#****************************************************************************
#*Table 1 and Supplementary Table 1 in manuscript
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
rm(list = ls())



message(paste0("Make sure the package \"e107\" \"gplots\" \"ROCR\" and \"multicore\" were installed!\n",
               "And \"Rgtsp\", which located in the AIMS package, were also installed!\n",
               "Besides, make sure there are no duplcated rownames in the expression matrix!!!"))


setwd("./Scripts/")
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

compare_array_RNA_Seq_pam_and_cor_pred <- function(training_and_tesing_res,expr_dat){
  model_used <- training_and_tesing_res
  #expr_dat <- expr_dat
  
  model_with_best_MCC <- sapply(model_used,function(x){
    mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
  }) %>% which.max()
  best_model <- model_used[[model_with_best_MCC]]
  testing_samples_in_RNA_seq <- intersect(colnames(expr_dat),best_model$predict_res$test_set %>% rownames())
  
  dat_for_pred <- expr_dat[best_model$training_model$centroids %>% rownames(),
                           testing_samples_in_RNA_seq]
  
  RNA_seq_pamr_pred <- pamr.predict(fit = best_model$training_model,
                                    newx = dat_for_pred,
                                    threshold = best_model$delt) %>% 
    as.data.frame() %>% mutate(barcode=testing_samples_in_RNA_seq) %>% rename(RNA_seq_pam_predict=1)
  
  RNA_seq_cor_pred <- cor(x=dat_for_pred[best_model$centroid %>% rownames(),],
                          y = best_model$centroid,method = "spearman") %>% 
    apply(1,function(x)names(which.max(x))) %>% 
    as.data.frame() %>% mutate(barcode=testing_samples_in_RNA_seq) %>% rename(RNA_seq_cor_predict=1)
  
  predicted_res <- reduce(list(best_model$predict_res$test_set %>% rownames_to_column("barcode"),
                               RNA_seq_pamr_pred,RNA_seq_cor_pred),inner_join,by="barcode")
  
  list(pam_compre=caret::confusionMatrix(data = factor(predicted_res$pam_predict,
                                                       levels = as.character(1:5)),
                                         reference=factor(predicted_res$RNA_seq_pam_predict,
                                                          levels = as.character(1:5))),
       cor_compare=caret::confusionMatrix(data = factor(predicted_res$cor_predict,
                                                        levels = as.character(1:5)),
                                          reference=factor(predicted_res$RNA_seq_cor_predict,
                                                           levels = as.character(1:5))),
       pred_res=predicted_res)
}


compare_array_RNA_Seq_ClaNC_PAM_Cor <- function(ClaNC_and_PAM_training_and_test_res,expr_dat){
  model_used <- ClaNC_and_PAM_training_and_test_res
  ClaNC_PAM_compre_res <- lapply(model_used,function(x)x$ClaNC_res)
  ClaNC_PAM_Cor_compare_res <- lapply(model_used,function(x)x$PAM_res)
  
  ClaNC_with_best_MCC <- sapply(ClaNC_PAM_compre_res,function(x){
    mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$table))
  }) %>% which.max()
  best_ClaNC_model <- ClaNC_PAM_compre_res[[ClaNC_with_best_MCC]]
  
  #The ClaNC_PAM and ClaNC_PAM_Cor_compare_res are using the same training and testing samples
  testing_samples_in_RNA_seq <- intersect(colnames(expr_dat),
                                          ClaNC_PAM_Cor_compare_res[[ClaNC_with_best_MCC]]$predict_res$test_set %>% 
                                            rownames())
  RNA_seq_ClaNC_predict <- predictClanc(data = expr_dat[,testing_samples_in_RNA_seq],geneNames = rownames(expr_dat),
                                        fit = best_ClaNC_model$predicting_model$builded)
  RNA_seq_ClaNC_predict <- best_ClaNC_model$class_label[RNA_seq_ClaNC_predict] %>% 
    setNames(testing_samples_in_RNA_seq) %>% as.data.frame() %>% rownames_to_column("barcode") %>% 
    rename(ClaNC_RNA_seq=2)
  
  
  array_ClaNC_predict <- setNames(best_ClaNC_model$testing_predicted_res,
                                  rownames(ClaNC_PAM_Cor_compare_res[[ClaNC_with_best_MCC]]$predict_res$test_set)) %>%
    as.data.frame(check.names=FALSE) %>% rownames_to_column("barcode") %>% rename(ClaNC_array=2)
  
  ClaNC_compare <- array_ClaNC_predict %>% inner_join(RNA_seq_ClaNC_predict,by = "barcode")
  
  Cla_PAM_and_Cla_PAM_Cor <- compare_array_RNA_Seq_pam_and_cor_pred(training_and_tesing_res=ClaNC_PAM_Cor_compare_res,
                                                                    expr_dat = expr_dat)
  list(ClaNC_compare=list(ClaNC_compare=caret::confusionMatrix(data = factor(ClaNC_compare$ClaNC_array,
                                                                             levels = as.character(1:5)),
                                                               reference=factor(ClaNC_compare$ClaNC_RNA_seq,
                                                                                levels = as.character(1:5))),
                          pred_res=ClaNC_compare),
       ClaNC_PAM_and_ClaNC_PAM_Cor_compare=Cla_PAM_and_Cla_PAM_Cor)
}


# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("../Data/")

load("./R_data/BRCA/Expr_and_pheno.RData")
expr_dat <- X1_genename
expr_dat <- log2(expr_dat+1)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]

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

RNA_seq_expr_for_clust <- tumor_expr


# array data
load("20221216_BRCA_microarray_expression.RData")

expr_dat <- BRCA_array_expr %>% apply(2,as.numeric)
rownames(expr_dat) <- rownames(BRCA_array_expr)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]

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

array_expr_for_clust <- tumor_expr


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


dat_for_model <- tumor_expr
res_prefix <- "TCGA_RNA_seq_and_array"
res_prefix_intrinsic <- "TCGA_RNA_seq_and_array_intrinsic"


setwd("../Results/")
#all intersect gene predicting model
load(paste0(res_prefix,"_AIMS_training_res.RData"))
a <- predict.one.vs.all.tsp(D = train_and_test_res$train_data,GeneName = rownames(train_and_test_res$train_data),
                            one.vs.all.tsp = train_and_test_res$training_res$final_model)

predicting_res <- a$cl %>% as.data.frame() %>% rownames_to_column("barcode") %>% rename(subtype=2) %>% 
  mutate(data_type=str_replace_all(barcode,"_TCGA.*","")) %>% 
  mutate(samples=str_replace_all(barcode,".*_","")) %>% 
  pivot_wider(id_cols = -barcode,names_from = data_type,values_from = subtype)

#**************************
#*Supplementary Table 1****
#*And its Kappa************
#**************************
table(predicting_res$RNA_seq,predicting_res$Array)
mltools::mcc(confusionM = as.matrix.data.frame(table(predicting_res$RNA_seq,predicting_res$Array)))
caret::confusionMatrix(data = factor(predicting_res$Array,levels = sort(unique(predicting_res$Array))),
                       reference=factor(predicting_res$RNA_seq,levels = sort(unique(predicting_res$RNA_seq))))

#********************************************** ****************************
# High consistency between RNA_seq and array may be attributed to  --------
# sample misusaging
#**************************************************************************




#************
#*Table 1****
#************

#predict RNA_seq data using AIMS by array data of all genes trained model
load("TCGA_array_intersect_gene_train_and_test_AIMS_training_res.RData")

AIMS_model <- train_and_test_res
expr_dat <- RNA_seq_expr

testing_samples_in_RNA_seq <- intersect(colnames(expr_dat),AIMS_model$predicted_res$cl %>% rownames())
array_testing_sam_class <- AIMS_model$predicted_res$cl[testing_samples_in_RNA_seq,] %>% 
  as.data.frame %>% rename(Cluster_in_array=1) %>% rownames_to_column("barcode")
RNA_seq_testing_sam_class <- predict.one.vs.all.tsp(D = expr_dat[,testing_samples_in_RNA_seq],
                                                    GeneName = rownames(expr_dat),
                                                    one.vs.all.tsp = AIMS_model$training_res$final_model)$cl %>% 
  as.data.frame %>% rename(Cluster_in_RNA_seq=1) %>% rownames_to_column("barcode")

combined_predicting <- array_testing_sam_class %>%  left_join(RNA_seq_testing_sam_class,by="barcode")

mltools::mcc(confusionM = as.matrix.data.frame(table(factor(combined_predicting$Cluster_in_array,
                                                            levels = as.character(1:5)),
                                                     factor(combined_predicting$Cluster_in_RNA_seq,
                                                            levels = as.character(1:5)))))

#********************
#*Gene-pairs all gene
#********************
caret::confusionMatrix(data = factor(combined_predicting$Cluster_in_array,levels = as.character(1:5)),
                       reference=factor(combined_predicting$Cluster_in_RNA_seq,levels = as.character(1:5)))



#predict RNA_seq data using AIMS by array data of intrinsic genes trained model
load("TCGA_array_intersect_gene_train_and_test_intrinsic_gene_AIMS_training_res.RData")

AIMS_model <- train_and_test_res_intrinsic
expr_dat <- RNA_seq_expr

testing_samples_in_RNA_seq <- intersect(colnames(expr_dat),AIMS_model$predicted_res$cl %>% rownames())
array_testing_sam_class <- AIMS_model$predicted_res$cl[testing_samples_in_RNA_seq,] %>% 
  as.data.frame %>% rename(Cluster_in_array=1) %>% rownames_to_column("barcode")
RNA_seq_testing_sam_class <- predict.one.vs.all.tsp(D = expr_dat[,testing_samples_in_RNA_seq],
                                                    GeneName = rownames(expr_dat),
                                                    one.vs.all.tsp = AIMS_model$training_res$final_model)$cl %>% 
  as.data.frame %>% rename(Cluster_in_RNA_seq=1) %>% rownames_to_column("barcode")

combined_predicting <- array_testing_sam_class %>%  left_join(RNA_seq_testing_sam_class,by="barcode")

mltools::mcc(confusionM = as.matrix.data.frame(table(factor(combined_predicting$Cluster_in_array,
                                                            levels = as.character(1:5)),
                                                     factor(combined_predicting$Cluster_in_RNA_seq,
                                                            levels = as.character(1:5)))))
#**************************
#*Gene-pairs intrinsic gene
#**************************
caret::confusionMatrix(data = factor(combined_predicting$Cluster_in_array,levels = as.character(1:5)),
                       reference=factor(combined_predicting$Cluster_in_RNA_seq,levels = as.character(1:5)))




#predict RNA_seq data using PAM by array data of all genes trained model
load("TCGA_array_intersect_gene_train_and_test_ClaNC_and_PAM_training_res.RData")
res <- compare_array_RNA_Seq_pam_and_cor_pred(training_and_tesing_res = PAM_and_PAM_plus_Cor,expr_dat = RNA_seq_expr)

#*************
#*PAM all gene
#*************
res$pam_compre


#**********************
#*PAM+Spearman all gene
#**********************
res$cor_compare


res <- compare_array_RNA_Seq_ClaNC_PAM_Cor(ClaNC_and_PAM_training_and_test_res = ClaNC_and_PAM_train,
                                           expr_dat = RNA_seq_expr)
#***************
#*ClaNC all gene
#***************
res$ClaNC_compare$ClaNC_compare

#*******************
#*ClaNC+PAM all gene
#*******************
res$ClaNC_PAM_and_ClaNC_PAM_Cor_compare$pam_compre

#****************************
#*ClaNC+PAM+Spearman all gene
#****************************
res$ClaNC_PAM_and_ClaNC_PAM_Cor_compare$cor_compare

#predict RNA_seq data using PAM by array data of intrinsic genes trained model
load("TCGA_array_intersect_gene_train_and_test_intrinsic_gene_ClaNC_and_PAM_training_res.RData")
res <- compare_array_RNA_Seq_pam_and_cor_pred(training_and_tesing_res = PAM_and_PAM_plus_Cor_intrinsic,
                                              expr_dat = RNA_seq_expr)

#*******************
#*PAM intrinsic gene
#*******************
res$pam_compre

#**********************
#*PAM+Spearman intrinsic gene
#**********************
res$cor_compare


res <- compare_array_RNA_Seq_ClaNC_PAM_Cor(ClaNC_and_PAM_training_and_test_res = ClaNC_and_PAM_train_intrinsic,
                                           expr_dat = RNA_seq_expr)
#***************
#*ClaNC intrinsic gene
#***************
res$ClaNC_compare$ClaNC_compare

#*******************
#*ClaNC+PAM intrinsic gene
#*******************
res$ClaNC_PAM_and_ClaNC_PAM_Cor_compare$pam_compre

#****************************
#*ClaNC+PAM+Spearman intrinsic gene
#****************************
res$ClaNC_PAM_and_ClaNC_PAM_Cor_compare$cor_compare


