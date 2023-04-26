##Chinese TNBC assignment
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
library(biomaRt)
#library(psych)


setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")



#*****************
#*Define functions
#*****************
#Cor_predict
cor_predict <- function(expr_dat,PAM_model){
  pamr_centroids <- pamr.predict(fit = PAM_model$training_model,newx = PAM_model$traning_and_testing$train_set$x,
                                 threshold = PAM_model$delt,type = "centroid")
  
  
  PAM_model$traning_and_testing$train_set$geneid <- 1:nrow(PAM_model$traning_and_testing$train_set$x)
  PAM_model$traning_and_testing$train_set$genenames <- rownames(PAM_model$traning_and_testing$train_set$x)
  pam_info <- capture.output(pamr_survived_genes <- pamr.listgenes(fit = PAM_model$training_model,
                                                                   data = PAM_model$traning_and_testing$train_set,
                                                                   threshold = PAM_model$delt,
                                                                   genenames = TRUE)[,"name"])
  if(length(pamr_survived_genes)<=1){
    message("Number of survived genes less or equal to 1 produced which would lead to error in the following analysis
            because structure of resulting data were changed!")
    return(NULL)
  }
  #pamr_centroids <- pamr_centroids[pamr_survived_genes,]
  
  pamr_centroids <- pamr_centroids[intersect(rownames(PAM_model$training_model$centroids),pamr_survived_genes),]
  
  pamr_centroids_test_cor <- cor(x = expr_dat[rownames(pamr_centroids),],y = pamr_centroids,method = "spearman")
  test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
}

predict_subtype <- function(tumor_expr_dat,training_and_testing_res,
                            predicting_type=c("AIMS","PAM","ClaNC"),class_label){
  switch(predicting_type,
         AIMS={
           res <- predict.one.vs.all.tsp(D = tumor_expr_dat,GeneName = rownames(tumor_expr_dat),
                                         one.vs.all.tsp = training_and_testing_res$training_res$final_model)$cl
           colnames(res)[1] <- "AIMS_predict"
           res %>% as.data.frame %>% rownames_to_column("barcode")
         },
         PAM={
           #Select PAM model with highest MCC
           predict_highest_MCC_PAM <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
           }))
           PAM_model <- training_and_testing_res[[predict_highest_MCC_PAM]]
           
           genes_in_model_not_in_data <- setdiff(rownames(PAM_model$training_model$centroids),rownames(tumor_expr_dat))
           if (length(genes_in_model_not_in_data) >=1) {
             message(paste0("Check if these genes existed in your data\n",
                            "or if these genes exist in your data with their alias!:\n",
                            paste0(genes_in_model_not_in_data,collapse = ", ")))
           }
           
           
           
           message(paste0("\nIf data contain NA data for some samples, you should fill these NA data!",
                          "\nif you not, *Error* might folow!"))
           
           all_NA_genes_in_data <- is.na(tumor_expr_dat[,1]) %>% which() %>% names()
           
           #drop NA genes the model which might affect the results
           Data_NA_genes_in_model <- intersect(all_NA_genes_in_data,rownames(PAM_model$training_model$centroids))
           
           #NA genes and genes in model not in data
           genes_NA_or_not_in_model <- c(genes_in_model_not_in_data,Data_NA_genes_in_model)
           
           if(length(genes_NA_or_not_in_model)>=1){
             PAM_model$training_model$centroids <- 
               PAM_model$training_model$centroids[setdiff(rownames(PAM_model$training_model$centroids),
                                                          genes_NA_or_not_in_model),]
             PAM_model$training_model$sd <- PAM_model$training_model$sd[setdiff(names(PAM_model$training_model$sd),
                                                                                genes_NA_or_not_in_model)]
             PAM_model$training_model$centroid.overall <- 
               PAM_model$training_model$centroid.overall[setdiff(names(PAM_model$training_model$centroid.overall),
                                                                 genes_NA_or_not_in_model)]
             
             message(paste0("\nThese genes were exclude from the PAM trained model:",
                            paste(genes_NA_or_not_in_model,collapse = ", "),
                            ",\nsince these genes were NA in expression data\n",
                            "THIS MIGHT CHANGE THE PREDICTING RESULTS!"))
           }
           
           
           PAM_predict <- pamr.predict(fit = PAM_model$training_model,
                                       newx = tumor_expr_dat[rownames(PAM_model$training_model$centroids),],
                                       threshold = PAM_model$delt)
           names(PAM_predict) <- colnames(tumor_expr_dat)
           
           #Select PAM + Cor model with highest MCC
           predict_highest_MCC_Cor <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
           }))
           PAM_model <- training_and_testing_res[[predict_highest_MCC_Cor]]
           
           genes_in_model_not_in_data <- setdiff(rownames(PAM_model$training_model$centroids),rownames(tumor_expr_dat))
           if (length(genes_in_model_not_in_data) >=1) {
             message(paste0("Check if these genes existed in your data\n",
                            "or if these genes exist in your data with their alias!:\n",
                            paste0(genes_in_model_not_in_data,collapse = ", ")))
           }
           
           
           
           message(paste0("\nIf data contain NA data for some samples, you should fill these NA data!",
                          "\nif you not, *Error* might folow!"))
           
           all_NA_genes_in_data <- is.na(tumor_expr_dat[,1]) %>% which() %>% names()
           
           #drop NA genes the model which might affect the results
           Data_NA_genes_in_model <- intersect(all_NA_genes_in_data,rownames(PAM_model$training_model$centroids))
           
           #NA genes and genes in model not in data
           genes_NA_or_not_in_model <- c(genes_in_model_not_in_data,Data_NA_genes_in_model)
           
           if(length(genes_NA_or_not_in_model)>=1){
             PAM_model$training_model$centroids <- 
               PAM_model$training_model$centroids[setdiff(rownames(PAM_model$training_model$centroids),
                                                          genes_NA_or_not_in_model),]
             PAM_model$training_model$sd <- PAM_model$training_model$sd[setdiff(names(PAM_model$training_model$sd),
                                                                                genes_NA_or_not_in_model)]
             PAM_model$training_model$centroid.overall <- 
               PAM_model$training_model$centroid.overall[setdiff(names(PAM_model$training_model$centroid.overall),
                                                                 genes_NA_or_not_in_model)]
             
             message(paste0("\nThese genes were exclude from the PAM trained model:",
                            paste(genes_NA_or_not_in_model,collapse = ", "),
                            ",\nsince these genes were NA in expression data\n",
                            "THIS MIGHT CHANGE THE PREDICTING RESULTS!"))
           }
           
           
           PAM_Cor_predict <- cor_predict(expr_dat = tumor_expr_dat,PAM_model = PAM_model)
           
           cbind(PAM_predict,PAM_Cor_predict) %>% as.data.frame %>% rownames_to_column("barcode")
         },
         ClaNC={
           if(missing(class_label))
             stop("Please providing class label when using ClaNC!")
           
           #Select ClaNC model with highest MCC
           predict_highest_MCC_ClaNC <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
           }))
           ClaNC_model <- training_and_testing_res[[predict_highest_MCC_ClaNC]]
           
           
           #genes in ClaNC model not in data
           genes_in_model_not_in_data_ClaNC <- setdiff(rownames(ClaNC_model$ClaNC_res$predicting_model$builded$cntrds),
                                                       rownames(tumor_expr_dat))
           
           
           if(length(genes_in_model_not_in_data_ClaNC) >= 1){
             message(paste0("\nThese genes were exclude from the ClaNC trained model:",
                            paste(genes_in_model_not_in_data_ClaNC,collapse = ", "),
                            ",\nsince these genes did not existed in expression data\n",
                            "THIS MIGHT CHANGE THE PREDICTING RESULTS!"))
             
             intersect_genes <- intersect(rownames(ClaNC_model$ClaNC_res$predicting_model$builded$cntrds),
                                          rownames(tumor_expr_dat))
             new_ClaNC_cntrds <- ClaNC_model$ClaNC_res$predicting_model$builded$cntrds[intersect_genes,]
             new_ClaNC_pooledSD <- ClaNC_model$ClaNC_res$predicting_model$builded$pooledSD[intersect_genes]
             ClaNC_model$ClaNC_res$predicting_model$builded$cntrds <- new_ClaNC_cntrds
             ClaNC_model$ClaNC_res$predicting_model$builded$geneNames <- intersect_genes
             ClaNC_model$ClaNC_res$predicting_model$builded$pooledSD <- new_ClaNC_pooledSD
           }
           
           
           
           clanc_res <- ClaNC_predicted(ClaNC_builded = ClaNC_model$ClaNC_res$predicting_model$builded,
                                        expr_set = list(x=tumor_expr_dat),
                                        gene_names = rownames(tumor_expr_dat),class_label = class_label)
           ClaNC_predict <- names(clanc_res)
           ClaNC_predict <- setNames(ClaNC_predict,colnames(tumor_expr_dat))
           
           #Select PAM model with highest MCC
           predict_highest_MCC_PAM <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
           }))
           PAM_model <- training_and_testing_res[[predict_highest_MCC_PAM]]$PAM_res
           
           genes_in_model_not_in_data <- setdiff(rownames(PAM_model$training_model$centroids),rownames(tumor_expr_dat))
           if (length(genes_in_model_not_in_data) >=1) {
             message(paste0("Check if these genes existed in your data\n",
                            "or if these genes exist in your data with their alias!:\n",
                            paste0(genes_in_model_not_in_data,collapse = ", ")))
           }
           
           
           
           message(paste0("\nIf data contain NA data for some samples, you should fill these NA data!",
                          "\nif you not, *Error* might folow!"))
           
           all_NA_genes_in_data <- is.na(tumor_expr_dat[,1]) %>% which() %>% names()
           
           #drop NA genes the model which might affect the results
           Data_NA_genes_in_model <- intersect(all_NA_genes_in_data,rownames(PAM_model$training_model$centroids))
           
           #NA genes and genes in model not in data
           genes_NA_or_not_in_model <- c(genes_in_model_not_in_data,Data_NA_genes_in_model)
           
           if(length(genes_NA_or_not_in_model)>=1){
             PAM_model$training_model$centroids <- 
               PAM_model$training_model$centroids[setdiff(rownames(PAM_model$training_model$centroids),
                                                          genes_NA_or_not_in_model),]
             PAM_model$training_model$sd <- PAM_model$training_model$sd[setdiff(names(PAM_model$training_model$sd),
                                                                                genes_NA_or_not_in_model)]
             PAM_model$training_model$centroid.overall <- 
               PAM_model$training_model$centroid.overall[setdiff(names(PAM_model$training_model$centroid.overall),
                                                                 genes_NA_or_not_in_model)]
             
             message(paste0("\nThese genes were exclude from the PAM trained model:",
                            paste(genes_NA_or_not_in_model,collapse = ", "),
                            ",\nsince these genes were NA in expression data\n",
                            "THIS MIGHT CHANGE THE PREDICTING RESULTS!"))
           }
           
           PAM_predict <- pamr.predict(fit = PAM_model$training_model,
                                       newx = tumor_expr_dat[rownames(PAM_model$training_model$centroids),],
                                       threshold = PAM_model$delt)
           names(PAM_predict) <- colnames(tumor_expr_dat)
           
           #Select PAM + Cor model with highest MCC
           predict_highest_MCC_Cor <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
           }))
           PAM_model <- training_and_testing_res[[predict_highest_MCC_Cor]]$PAM_res
           
           genes_in_model_not_in_data <- setdiff(rownames(PAM_model$training_model$centroids),rownames(tumor_expr_dat))
           if (length(genes_in_model_not_in_data) >=1) {
             message(paste0("Check if these genes existed in your data\n",
                            "or if these genes exist in your data with their alias!:\n",
                            paste0(genes_in_model_not_in_data,collapse = ", ")))
           }
           
           
           
           message(paste0("\nIf data contain NA data for some samples, you should fill these NA data!",
                          "\nif you not, *Error* might folow!"))
           
           all_NA_genes_in_data <- is.na(tumor_expr_dat[,1]) %>% which() %>% names()
           
           #drop NA genes the model which might affect the results
           Data_NA_genes_in_model <- intersect(all_NA_genes_in_data,rownames(PAM_model$training_model$centroids))
           
           #NA genes and genes in model not in data
           genes_NA_or_not_in_model <- c(genes_in_model_not_in_data,Data_NA_genes_in_model)
           
           if(length(genes_NA_or_not_in_model)>=1){
             PAM_model$training_model$centroids <- 
               PAM_model$training_model$centroids[setdiff(rownames(PAM_model$training_model$centroids),
                                                          genes_NA_or_not_in_model),]
             PAM_model$training_model$sd <- PAM_model$training_model$sd[setdiff(names(PAM_model$training_model$sd),
                                                                                genes_NA_or_not_in_model)]
             PAM_model$training_model$centroid.overall <- 
               PAM_model$training_model$centroid.overall[setdiff(names(PAM_model$training_model$centroid.overall),
                                                                 genes_NA_or_not_in_model)]
             
             message(paste0("\nThese genes were exclude from the PAM trained model:",
                            paste(genes_NA_or_not_in_model,collapse = ", "),
                            ",\nsince these genes were NA in expression data\n",
                            "THIS MIGHT CHANGE THE PREDICTING RESULTS!"))
           }
           
           
           PAM_Cor_predict <- cor_predict(expr_dat = tumor_expr_dat,PAM_model = PAM_model)
           
           cbind(ClaNC_predict,ClaNC_PAM_predict=PAM_predict,ClaNC_PAM_Cor_predict=PAM_Cor_predict) %>% 
             as.data.frame %>% rownames_to_column("barcode")
         })
  
}



setwd("/mnt/Miscrosoft/TCGA_project/Data/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/")
setwd("./Gene_Expression_Quantification/0019c951-16c5-48d0-85c8-58d96b12d330/")
Ensembl_and_Entrez <- read.table("ba295155-272e-43eb-9d6a-e4c9c392e68b.rna_seq.augmented_star_gene_counts.tsv",
                                 header = TRUE,
                                 sep = "\t",stringsAsFactors = FALSE)[-1:-4,1:3]
Ensembl_and_Entrez$gene_id <- str_remove_all(Ensembl_and_Entrez$gene_id,"\\..*")

setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/AIMS_assumption/")
Chinese_TNBC_rep_samples_new <- read.table("Data/TNBC_43samples_genes_fpkm.csv",
                                           header = TRUE,sep = ",",stringsAsFactors = FALSE)
rownames(Chinese_TNBC_rep_samples_new) <- Chinese_TNBC_rep_samples_new[,1]
Chinese_TNBC_rep_samples_new <- Chinese_TNBC_rep_samples_new[,-1]
# data_to_PAM50_new <- log2(Chinese_TNBC_rep_samples_new+0.01)



rep_samples <- colnames(Chinese_TNBC_rep_samples_new)
C_TNBC_expr <- apply(Chinese_TNBC_rep_samples_new[Ensembl_and_Entrez[,1],rep_samples],
                      2,function(x,y){
                        tapply(x,y,mean,na.rm=TRUE)
                      },y=Ensembl_and_Entrez[,2])

C_TNBC_expr <- log2(C_TNBC_expr+0.01)


setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
#RNA-seq All intersect genes predicting res
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
load("../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")
load("TCGA_BRCA_clinical_dat.RData")


#AIMS
all_gene_AIMS_res <- predict_subtype(tumor_expr_dat = C_TNBC_expr,
                                     training_and_testing_res = train_and_test_res,
                                     predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res <- predict_subtype(tumor_expr_dat = C_TNBC_expr,training_and_testing_res = PAM_and_PAM_plus_Cor,
                                    predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = C_TNBC_expr,training_and_testing_res = ClaNC_and_PAM_train,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))






setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/AIMS_assumption/")
Chinese_TNBC_rep_samples_old <- read.table("Data/tnbc2019_448samples_dir_fpkm.csv",header = TRUE,
                                           sep = ",",stringsAsFactors = FALSE)
rownames(Chinese_TNBC_rep_samples_old) <- Chinese_TNBC_rep_samples_old[,1]
Chinese_TNBC_rep_samples_old <- Chinese_TNBC_rep_samples_old[,-1]


C_TNBC_expr_old <- apply(Chinese_TNBC_rep_samples_old[Ensembl_and_Entrez[,1],rep_samples],
                         2,function(x,y){
                           tapply(x,y,mean,na.rm=TRUE)
                         },y=Ensembl_and_Entrez[,2])

C_TNBC_expr_old <- log2(C_TNBC_expr_old+0.01)


#AIMS
all_gene_AIMS_res_old <- predict_subtype(tumor_expr_dat = C_TNBC_expr_old,
                                         training_and_testing_res = train_and_test_res,
                                         predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res_old <- predict_subtype(tumor_expr_dat = C_TNBC_expr_old,
                                        training_and_testing_res = PAM_and_PAM_plus_Cor,
                                        predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res_old <- predict_subtype(tumor_expr_dat = C_TNBC_expr_old,
                                          training_and_testing_res = ClaNC_and_PAM_train,
                                          predicting_type = "ClaNC",
                                          class_label = levels(as.factor(datset_for_model$y)))


#All genes predict results
all_genes_predict_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res,
                                     all_gene_AIMS_res_old,all_gene_PAM_res_old,all_gene_ClaNC_res_old),
                                left_join,by="barcode",suffix=c(".new",".old"))

all_genes_kappa <- sapply(str_remove_all(colnames(all_genes_predict_res[2:7]),"\\.new"),function(x,y){
  caret::confusionMatrix(data = factor(all_genes_predict_res[,paste0(x,".new")],levels = y),
                         reference=factor(all_genes_predict_res[,paste0(x,".old")],levels = y))$overall["Kappa"]
},y=as.character(1:5),simplify = TRUE)


# all_genes_kappa <- sapply(str_remove_all(colnames(all_genes_predict_res[2:7]),"\\.new"),function(x,y){
#   psych::cohen.kappa(x=all_genes_predict_res[,paste0(x,c(".new",".old"))])
# },y=as.character(1:5),simplify = FALSE)




setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
#RNA-seq Intrinsic intersect genes predicting res
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData")
load("../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")
load("TCGA_BRCA_clinical_dat.RData")


#AIMS
Intrinsic_gene_AIMS_res <- predict_subtype(tumor_expr_dat = C_TNBC_expr,
                                           training_and_testing_res = train_and_test_res_intrinsic,
                                           predicting_type = "AIMS")
#PAM and PAM + Cor
Intrinsic_gene_PAM_res <- predict_subtype(tumor_expr_dat = C_TNBC_expr,
                                          training_and_testing_res = PAM_and_PAM_plus_Cor_intrinsic,
                                          predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
Intrinsic_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = C_TNBC_expr,
                                            training_and_testing_res = ClaNC_and_PAM_train_intrinsic,
                                            predicting_type = "ClaNC",
                                            class_label = levels(as.factor(datset_for_model$y)))


#AIMS
Intrinsic_gene_AIMS_res_old <- predict_subtype(tumor_expr_dat = C_TNBC_expr_old,
                                               training_and_testing_res = train_and_test_res_intrinsic,
                                               predicting_type = "AIMS")
#PAM and PAM + Cor
Intrinsic_gene_PAM_res_old <- predict_subtype(tumor_expr_dat = C_TNBC_expr_old,
                                              training_and_testing_res = PAM_and_PAM_plus_Cor_intrinsic,
                                              predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
Intrinsic_gene_ClaNC_res_old <- predict_subtype(tumor_expr_dat = C_TNBC_expr_old,
                                                training_and_testing_res = ClaNC_and_PAM_train_intrinsic,
                                                predicting_type = "ClaNC",
                                                class_label = levels(as.factor(datset_for_model$y)))



Intrinsic_genes_predict_res <- reduce(list(Intrinsic_gene_AIMS_res,Intrinsic_gene_PAM_res,Intrinsic_gene_ClaNC_res,
                                           Intrinsic_gene_AIMS_res_old,Intrinsic_gene_PAM_res_old,
                                           Intrinsic_gene_ClaNC_res_old),
                                      left_join,by="barcode",suffix=c(".new",".old"))



Intrinsic_genes_kappa <- sapply(str_remove_all(colnames(Intrinsic_genes_predict_res[2:7]),"\\.new"),function(x,y){
  caret::confusionMatrix(data = factor(Intrinsic_genes_predict_res[,paste0(x,".new")],levels = y),
                         reference=factor(Intrinsic_genes_predict_res[,paste0(x,".old")],levels = y))$overall["Kappa"]
},y=as.character(1:5),simplify = TRUE)


# Intrinsic_genes_kappa <- sapply(str_remove_all(colnames(Intrinsic_genes_predict_res[2:7]),"\\.new"),function(x,y){
#   psych::cohen.kappa(x=Intrinsic_genes_predict_res[,paste0(x,c(".new",".old"))])
# },y=as.character(1:5),simplify = FALSE)


all_intrinsic_predict_res <- left_join(all_genes_predict_res,Intrinsic_genes_predict_res,
                                       by="barcode",suffix=c(".all",".intrinsic"))



all_Intrinsic_genes_kappa <- sapply(str_remove_all(colnames(all_intrinsic_predict_res[2:7]),"\\.new.*"),function(x,y){
  caret::confusionMatrix(data = factor(all_intrinsic_predict_res[,paste0(x,".new.all")],levels = y),
                         reference=factor(all_intrinsic_predict_res[,paste0(x,".new.intrinsic")],
                                          levels = y))$overall["Kappa"]
},y=as.character(1:5),simplify = TRUE)


all_Intrinsic_genes_kappa_old <- sapply(str_remove_all(colnames(all_intrinsic_predict_res[2:7]),"\\.new.*"),
                                        function(x,y){
  caret::confusionMatrix(data = factor(all_intrinsic_predict_res[,paste0(x,".old.all")],levels = y),
                         reference=factor(all_intrinsic_predict_res[,paste0(x,".old.intrinsic")],
                                          levels = y))$overall["Kappa"]
},y=as.character(1:5),simplify = TRUE)






setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
# intrinsic_genes <- read.table("2009_JCO_intrinsic_genes_S_table_5.txt",header = FALSE,sep = "\t",
#                               stringsAsFactors = FALSE)
# 
# 
# C_TNBC_dataset <- list(expr=cbind(C_TNBC_expr[intersect(rownames(C_TNBC_expr),intrinsic_genes[,1]),],
#                                   C_TNBC_expr_old[intersect(rownames(C_TNBC_expr_old),intrinsic_genes[,1]),]),
#                        batch=rep(c("batch_1","batch_2"),each=43))

C_TNBC_dataset <- list(expr=cbind(C_TNBC_expr,C_TNBC_expr_old),batch=rep(c("batch_1","batch_2"),each=43))
colnames(C_TNBC_dataset$expr) <- paste0(colnames(C_TNBC_dataset$expr),".",rep(c("batch_1","batch_2"),each=43))


C_TNBC_dataset$expr <- C_TNBC_dataset$expr[(apply(C_TNBC_dataset$expr[,1:43],1,sd,na.rm=TRUE)!=0)&
                                             (apply(C_TNBC_dataset$expr[,-(1:43)],1,sd,na.rm=TRUE)!=0),]
C_TNBC_dataset$expr <- C_TNBC_dataset$expr[!rownames(C_TNBC_dataset$expr) %>% is.na,]



combined_dat <- C_TNBC_dataset$expr
# combined_dat <- t(apply(combined_dat,1,function(x){
#   (x-mean(x))/sd(x)
# }))
annotation_col <- data.frame(Batch=C_TNBC_dataset$batch,row.names = colnames(combined_dat))

pdf("../Results/Figure/Figure_5_supp_original_value.pdf",height = 10,width = 10)
#pheatmap(mat = combined_dat[sample(1:nrow(combined_dat),size = 10000,replace = FALSE),],
pheatmap(mat = combined_dat,
         labels_row = rep('',nrow(combined_dat)),
         labels_col = colnames(combined_dat),
         annotation_col = annotation_col,fontsize_col = 2,
         clustering_method = "ward.D",
         #annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))
dev.off()

# PCA <- prcomp(x = t(combined_dat))
# PCA$rotation %>% dim()
# plot(PCA$rotation[,1],PCA$rotation[,2],col=rep(c("red","blue"),each=43))
# 
# PCA.sum <- summary(PCA)
# plot(PCA.sum$x[,1:2],col=rep(c("red","blue"),each=43))




#Z-Sccore separately
a <- t(apply(C_TNBC_dataset$expr[,1:43],1,function(x)(x-mean(x))/sd(x)))
b <- t(apply(C_TNBC_dataset$expr[,-(1:43)],1,function(x)(x-mean(x))/sd(x)))

combined_dat <- cbind(a,b)

# combined_dat <- t(apply(combined_dat,1,function(x){
#   (x-mean(x))/sd(x)
# }))
annotation_col <- data.frame(Batch=C_TNBC_dataset$batch,row.names = colnames(combined_dat))

pdf("../Results/Figure/Figure_5_supp_scaled_separatly.pdf",height = 10,width = 10)
# pheatmap(mat = combined_dat[sample(1:nrow(combined_dat),size = 10000,replace = FALSE),],
pheatmap(mat = combined_dat,
         labels_row = rep('',nrow(combined_dat)),
         labels_col = colnames(combined_dat),
         annotation_col = annotation_col,fontsize_col = 1,clustering_method = "ward.D",
         #annotation_colors = ann_colors,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))
dev.off()





# svg("../Results/Figure/Figure_5_supp_scaled_separatly.svg",height = 10,width = 10)
# pheatmap(mat = combined_dat[sample(1:nrow(combined_dat),size = 10000,replace = FALSE),],
# # pheatmap(mat = combined_dat,
#          labels_row = rep('',nrow(combined_dat)),
#          labels_col = colnames(combined_dat),
#          annotation_col = annotation_col,fontsize_col = 1,clustering_method = "ward.D",
#          #annotation_colors = ann_colors,
#          color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))
# dev.off()
# 





