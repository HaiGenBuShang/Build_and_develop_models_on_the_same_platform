library(TCGAbiolinks)
library(SummarizedExperiment)
library(sigclust2)
library(tidyverse)
library(dendextend)
library(pheatmap)
library(pamr)
library(survminer)
library(seqc)
rm(list = ls())

setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")


#*****************
#*Define functions
#*****************
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}



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


obtain_loaded_dat <- function(RData_file){
  load(RData_file)
  rm(RData_file)
  obj <- ls()
  #rm(RData_file)
  sapply(obj,function(x)list(get(x)))
}

#RData_file <- "../Results/TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData"

kappa_between_array_and_RNA_seq <- function(RNA_pred,array_pred,
                                            compare_part="$train_and_test_res$predicted_res$cl[,1]"){
  RNA_seq_subtype <- factor(eval(parse(text = paste0("RNA_pred",compare_part))),levels = as.character(1:5))
  array_subtype <- factor(eval(parse(text = paste0("array_pred",compare_part))),levels = as.character(1:5))
  caret::confusionMatrix(data = RNA_seq_subtype,reference=array_subtype)
}

# RNA_seq_intersect <- obtain_loaded_dat(paste0("../Results/",
#                                               "TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData"))
# array_intersect <- obtain_loaded_dat("../Results/TCGA_array_intersect_gene_consensus_sample_AIMS_training_res.RData")
# kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect)




#using PAM50 
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")



#***********************************
#*RNA-seq data model predict RNA-seq
#***********************************
#AIMS all intersect genes
load("../Results/SEQC_AGR_BGI_Z_score_AGR_model_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))




#all intersect genes
load("../Results/SEQC_AGR_BGI_Z_score_AGR_model_ClaNC_and_PAM_training_res.RData")
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



#RNA-seq predicting results for test sample of each predicting strategy










#*****************************
#**********predict consistency
#*****************************

setwd("/mnt/Miscrosoft/TCGA_project/Data/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/")
setwd("./Gene_Expression_Quantification/0019c951-16c5-48d0-85c8-58d96b12d330/")
Ensembl_and_Entrez <- read.table("ba295155-272e-43eb-9d6a-e4c9c392e68b.rna_seq.augmented_star_gene_counts.tsv",
                                 header = TRUE,
                                 sep = "\t",stringsAsFactors = FALSE)[-1:-4,1:3]
Ensembl_and_Entrez$gene_id <- str_remove_all(Ensembl_and_Entrez$gene_id,"\\..*")


SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_BGI[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_BGI[,3])


SEQC_ILM_BGI_FPKM_symbol <- apply(SEQC_ILM_BGI_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_BGI[,2])


SEQC_expr_BGI <- log2(SEQC_ILM_BGI_FPKM_symbol[intersect(Ensembl_and_Entrez$gene_name,
                                                         rownames(SEQC_ILM_BGI_FPKM_symbol)),]+0.01)


SEQC_ILM_CNL_FPKM_entrez <- apply(ILM_refseq_gene_CNL[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_CNL[,3])


SEQC_ILM_CNL_FPKM_symbol <- apply(SEQC_ILM_CNL_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_CNL[,2])


SEQC_expr_CNL <- log2(SEQC_ILM_CNL_FPKM_symbol[intersect(Ensembl_and_Entrez$gene_name,
                                                         rownames(SEQC_ILM_CNL_FPKM_symbol)),]+0.01)


setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
load("../Results/SEQC_AGR_BGI_Z_score_AGR_model_AIMS_training_res.RData")
load("../Results/SEQC_AGR_BGI_Z_score_AGR_model_ClaNC_and_PAM_training_res.RData")
load("../Results/SEQC_datset_for_model_and_dataset_for_model_intrinsic_g.RData")

#AIMS
all_gene_AIMS_res <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                     training_and_testing_res = train_and_test_res,
                                     predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,training_and_testing_res = PAM_and_PAM_plus_Cor,
                                    predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,training_and_testing_res = ClaNC_and_PAM_train,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))



#AIMS
all_gene_AIMS_res_CNL <- predict_subtype(tumor_expr_dat = SEQC_expr_CNL,
                                         training_and_testing_res = train_and_test_res,
                                         predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res_CNL <- predict_subtype(tumor_expr_dat = SEQC_expr_CNL,
                                        training_and_testing_res = PAM_and_PAM_plus_Cor,
                                        predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res_CNL <- predict_subtype(tumor_expr_dat = SEQC_expr_CNL,
                                          training_and_testing_res = ClaNC_and_PAM_train,
                                          predicting_type = "ClaNC",
                                          class_label = levels(as.factor(datset_for_model$y)))


#All genes predict results
all_genes_predict_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res,
                                     all_gene_AIMS_res_CNL,all_gene_PAM_res_CNL,all_gene_ClaNC_res_CNL),
                                inner_join,by="barcode",suffix=c(".BGI",".CNL")) %>% 
  filter(!grepl("^E|^F",barcode))



all_genes_kappa <- sapply(str_remove_all(colnames(all_genes_predict_res[2:7]),"\\.BGI"),function(x,y){
  caret::confusionMatrix(data = factor(all_genes_predict_res[,paste0(x,".BGI")],levels = y),
                         reference=factor(all_genes_predict_res[,paste0(x,".CNL")],levels = y))$overall["Kappa"]
},y=as.character(1:4),simplify = TRUE)




load("../Results/SEQC_AGR_BGI_Z_score_AGR_model_high_expr_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))













b <- PAM_and_PAM_plus_Cor$`The 1st errors`$traning_and_testing$train_set
b$geneid <- 1:nrow(b$x)
b$genenames <- rownames(b$x)
a <- pamr.listgenes(fit = PAM_and_PAM_plus_Cor$`The 1st errors`$training_model,
                    data = b,threshold = 57.41326,genenames = TRUE)[,"name"]
b <- pamr.predict(fit = PAM_and_PAM_plus_Cor$`The 1st errors`$training_model,threshold = 57.41326,type = "centroid")

b[a,]



MCC_for_BGI_and_CNL <- all_genes_predict_res

MCC_for_BGI_and_CNL[,-1] <- apply(MCC_for_BGI_and_CNL[,-1],2,function(x){
  x %>% factor(levels = as.character(1:4)) %>% as.integer()
})
MCC_for_BGI_and_CNL <- MCC_for_BGI_and_CNL %>% mutate(real_type = str_sub(barcode,1,1)) %>% 
  mutate(real_type=factor(real_type,levels = LETTERS[1:4]) %>% as.integer())


# MCC_for_BGI_and_CNL <- MCC_for_BGI_and_CNL %>% 
#   pivot_longer(cols = c(2,8),names_to = "AIMS_type",values_to = "AIMS") %>% 
#   pivot_longer(cols = c(2,7),names_to = "PAM_type",values_to = "PAM") %>% 
#   pivot_longer(cols = c(2,6),names_to = "PAM_cor_type",values_to = "PAM_Cor") %>% 
#   pivot_longer(cols = c(2,5),names_to = "ClaNC_type",values_to = "ClaNC") %>% 
#   pivot_longer(cols = c(2,4),names_to = "ClaNC_PAM_type",values_to = "ClaNC_PAM") %>% 
#   pivot_longer(cols = c(2,3),names_to = "ClaNC_PAM_Cor_type",values_to = "ClaNC_PAM_Cor")

MCC_for_BGI_and_CNL <- MCC_for_BGI_and_CNL %>%
  pivot_longer(cols = c(2:13))


sapply(str_remove_all(colnames(all_genes_predict_res[2:7]),"\\.BGI"),function(x){
  dat <- MCC_for_BGI_and_CNL %>% dplyr::filter(grepl(paste0("^",x),name))
  mltools::mcc(preds = dat$value,actuals=dat$real_type)
})


all_genes_kappa <- sapply(str_remove_all(colnames(all_genes_predict_res[2:7]),"\\.BGI"),function(x,y){
  caret::confusionMatrix(data = factor(all_genes_predict_res[,paste0(x,".BGI")],levels = y),
                         reference=factor(all_genes_predict_res[,paste0(x,".CNL")],levels = y))$overall["Kappa"]
},y=as.character(1:4),simplify = TRUE)



#cor predict using Pearson correlation
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
  
  pamr_centroids_test_cor <- cor(x = expr_dat[rownames(pamr_centroids),],y = pamr_centroids)
  test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
}



#AIMS
all_gene_AIMS_res <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                     training_and_testing_res = train_and_test_res,
                                     predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,training_and_testing_res = PAM_and_PAM_plus_Cor,
                                    predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,training_and_testing_res = ClaNC_and_PAM_train,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))



#AIMS
all_gene_AIMS_res_CNL <- predict_subtype(tumor_expr_dat = SEQC_expr_CNL,
                                         training_and_testing_res = train_and_test_res,
                                         predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res_CNL <- predict_subtype(tumor_expr_dat = SEQC_expr_CNL,
                                        training_and_testing_res = PAM_and_PAM_plus_Cor,
                                        predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res_CNL <- predict_subtype(tumor_expr_dat = SEQC_expr_CNL,
                                          training_and_testing_res = ClaNC_and_PAM_train,
                                          predicting_type = "ClaNC",
                                          class_label = levels(as.factor(datset_for_model$y)))


#All genes predict results
all_genes_predict_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res,
                                     all_gene_AIMS_res_CNL,all_gene_PAM_res_CNL,all_gene_ClaNC_res_CNL),
                                inner_join,by="barcode",suffix=c(".BGI",".CNL")) %>% 
  filter(!grepl("^E|^F",barcode))



MCC_for_BGI_and_CNL <- all_genes_predict_res

MCC_for_BGI_and_CNL[,-1] <- apply(MCC_for_BGI_and_CNL[,-1],2,function(x){
  x %>% factor(levels = as.character(1:4)) %>% as.integer()
})
MCC_for_BGI_and_CNL <- MCC_for_BGI_and_CNL %>% mutate(real_type = str_sub(barcode,1,1)) %>% 
  mutate(real_type=factor(real_type,levels = LETTERS[1:4]) %>% as.integer())


MCC_for_BGI_and_CNL <- MCC_for_BGI_and_CNL %>%
  pivot_longer(cols = c(2:13))


sapply(str_remove_all(colnames(all_genes_predict_res[2:7]),"\\.BGI"),function(x){
  dat <- MCC_for_BGI_and_CNL %>% dplyr::filter(grepl(paste0("^",x),name))
  mltools::mcc(preds = dat$value,actuals=dat$real_type)
})
