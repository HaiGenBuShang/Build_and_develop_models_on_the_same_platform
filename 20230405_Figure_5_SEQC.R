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
library(seqc)
library(dendextend)



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


SEQC_ILM_AGR_FPKM_entrez <- apply(ILM_refseq_gene_AGR[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_AGR[,3])

# SEQC_ILM_AGR_FPKM_entrez <- apply(ILM_refseq_gene_COH[,-1:-4],2,function(x,y){
#   x/((y/1000)*(sum(x)/10e5))
# },y=ILM_refseq_gene_COH[,3])





SEQC_ILM_AGR_FPKM_symbol <- apply(SEQC_ILM_AGR_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_AGR[,2])

SEQC_expr_AGR <- log2(SEQC_ILM_AGR_FPKM_symbol[intersect(Ensembl_and_Entrez$gene_name,
                                                         rownames(SEQC_ILM_AGR_FPKM_symbol)),]+0.01)




setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
#RNA-seq All intersect genes predicting res
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
load("../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")



#AIMS
all_gene_AIMS_res <- predict_subtype(tumor_expr_dat = SEQC_expr_AGR,
                                     training_and_testing_res = train_and_test_res,
                                     predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res <- predict_subtype(tumor_expr_dat = SEQC_expr_AGR,training_and_testing_res = PAM_and_PAM_plus_Cor,
                                    predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = SEQC_expr_AGR,training_and_testing_res = ClaNC_and_PAM_train,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))





SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_BGI[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_BGI[,3])

# SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_CNL[,-1:-4],2,function(x,y){
#   x/((y/1000)*(sum(x)/10e5))
# },y=ILM_refseq_gene_CNL[,3])



SEQC_ILM_BGI_FPKM_symbol <- apply(SEQC_ILM_BGI_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_BGI[,2])


SEQC_expr_BGI <- log2(SEQC_ILM_BGI_FPKM_symbol[intersect(Ensembl_and_Entrez$gene_name,
                                                         rownames(SEQC_ILM_BGI_FPKM_symbol)),]+0.01)


#AIMS
all_gene_AIMS_res_old <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                         training_and_testing_res = train_and_test_res,
                                         predicting_type = "AIMS")
#PAM and PAM + Cor
all_gene_PAM_res_old <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                        training_and_testing_res = PAM_and_PAM_plus_Cor,
                                        predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res_old <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                          training_and_testing_res = ClaNC_and_PAM_train,
                                          predicting_type = "ClaNC",
                                          class_label = levels(as.factor(datset_for_model$y)))


#All genes predict results
all_genes_predict_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res,
                                     all_gene_AIMS_res_old,all_gene_PAM_res_old,all_gene_ClaNC_res_old),
                                inner_join,by="barcode",suffix=c(".new",".old")) #%>% 
  #filter(grepl("^A",barcode))




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
#load("TCGA_BRCA_clinical_dat.RData")


#AIMS
Intrinsic_gene_AIMS_res <- predict_subtype(tumor_expr_dat = SEQC_expr_AGR,
                                           training_and_testing_res = train_and_test_res_intrinsic,
                                           predicting_type = "AIMS")
#PAM and PAM + Cor
Intrinsic_gene_PAM_res <- predict_subtype(tumor_expr_dat = SEQC_expr_AGR,
                                          training_and_testing_res = PAM_and_PAM_plus_Cor_intrinsic,
                                          predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
Intrinsic_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = SEQC_expr_AGR,
                                            training_and_testing_res = ClaNC_and_PAM_train_intrinsic,
                                            predicting_type = "ClaNC",
                                            class_label = levels(as.factor(datset_for_model$y)))


#AIMS
Intrinsic_gene_AIMS_res_old <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                               training_and_testing_res = train_and_test_res_intrinsic,
                                               predicting_type = "AIMS")
#PAM and PAM + Cor
Intrinsic_gene_PAM_res_old <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                              training_and_testing_res = PAM_and_PAM_plus_Cor_intrinsic,
                                              predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
Intrinsic_gene_ClaNC_res_old <- predict_subtype(tumor_expr_dat = SEQC_expr_BGI,
                                                training_and_testing_res = ClaNC_and_PAM_train_intrinsic,
                                                predicting_type = "ClaNC",
                                                class_label = levels(as.factor(datset_for_model$y)))



Intrinsic_genes_predict_res <- reduce(list(Intrinsic_gene_AIMS_res,Intrinsic_gene_PAM_res,Intrinsic_gene_ClaNC_res,
                                           Intrinsic_gene_AIMS_res_old,Intrinsic_gene_PAM_res_old,
                                           Intrinsic_gene_ClaNC_res_old),
                                      inner_join,by="barcode",suffix=c(".new",".old"))



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
# C_TNBC_dataset <- list(expr=cbind(SEQC_expr_AGR[intersect(rownames(SEQC_expr_AGR),intrinsic_genes[,1]),],
#                                   SEQC_expr_BGI[intersect(rownames(SEQC_expr_BGI),intrinsic_genes[,1]),]),
#                        batch=rep(c("AGR","BGI"),times=c(ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI))))

C_TNBC_dataset <- list(expr=cbind(SEQC_expr_AGR,SEQC_expr_BGI),
                       batch=rep(c("AGR","BGI"),times=c(ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI))))
colnames(C_TNBC_dataset$expr) <- paste0(colnames(C_TNBC_dataset$expr),".",
                                        rep(c("AGR","BGI"),times=c(ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI))))


# C_TNBC_dataset$expr <- C_TNBC_dataset$expr[(apply(C_TNBC_dataset$expr[,1:43],1,sd,na.rm=TRUE)!=0)&
#                                              (apply(C_TNBC_dataset$expr[,-(1:43)],1,sd,na.rm=TRUE)!=0),]
# C_TNBC_dataset$expr <- C_TNBC_dataset$expr[!rownames(C_TNBC_dataset$expr) %>% is.na,]



combined_dat <- C_TNBC_dataset$expr
# combined_dat <- t(apply(combined_dat,1,function(x){
#   (x-mean(x))/sd(x)
# }))
annotation_col <- data.frame(Batch=C_TNBC_dataset$batch,row.names = colnames(combined_dat),
                             Samples=colnames(C_TNBC_dataset$expr) %>% str_sub(1,1))


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


clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D") %>% as.dendrogram()




clustering_res <- hclust(d = dist(t(combined_dat[,grepl("^A",colnames(combined_dat))])),
                         method = "ward.D") %>% as.dendrogram()
clustering_res %>%  
  # dendextend::set("branches_col",
  #                 ifelse(grepl("AGR",labels(clustering_res)),"#E41A1C","#377EB8")) %>%
  set("labels_colors",annotation_col[labels(clustering_res),"Samples"] %>% 
        factor(levels = LETTERS[1:6]) %>% as.integer()) %>% 
  set("labels_cex",1) %>% 
  # set("branches_k_color",1) %>% 
  assign_values_to_leaves_edgePar(value = ifelse(grepl("AGR",labels(clustering_res)),"#E41A1C","#377EB8"),
                                  edgePar = "col") %>% 
  plot(main="Unnormalized Sample A data",cex.main=2)
  
legend("topright",legend = c("AGR","BGI"),
       col=c("#E41A1C","#377EB8"),
             lwd=1,bty = "n",cex = 1.5)








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






##********************************************************
##*Assess why samples were predicted as different subtypes
##********************************************************


all_genes_predict_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res,
                                     all_gene_AIMS_res_old,all_gene_PAM_res_old,all_gene_ClaNC_res_old),
                                inner_join,by="barcode",suffix=c(".new",".old"))



all_genes_predict_res[,c(1,2,8)]

train_and_test_res$training_res$final_model$all.pairs
"HHAT<UXS1" 
"NARF<PPOX"


train_and_test_res$training_res$final_model$selected.pairs.list$`3`
"OSBPL1A<ORMDL2"


all_genes_predict_res[all_genes_predict_res$AIMS_predict.new!=all_genes_predict_res$AIMS_predict.old,c(1,2,8)]

cbind(SEQC_expr_AGR[c("OSBPL1A","ORMDL2"),"A_1_L01_FlowCellA"],SEQC_expr_BGI[c("OSBPL1A","ORMDL2"),"A_1_L01_FlowCellA"])




find_expr_direction <- function(gene_order="OSBPL1A<ORMDL2",sample_expr){
  genes <- str_split_fixed(gene_order,pattern = "<",n=2)[1,]
  direction <- sample_expr[genes[1]]>sample_expr[genes[2]]
  names(direction) <- NULL
  ifelse(direction,paste0(genes,collapse = ">"),paste0(genes,collapse = "<"))
}


compare_two_sample <- function(gene_order,sample_expr_1,sample_expr_2){
  direction_sample_1 <- find_expr_direction(gene_order = gene_order,sample_expr = sample_expr_1)
  direction_sample_2 <- find_expr_direction(gene_order = gene_order,sample_expr = sample_expr_2)
  direction_sample_1==direction_sample_2
}


compare_two_sample(gene_order = "ITGB8<NEIL3",
                   sample_expr_1 = SEQC_expr_AGR[,"A_1_L01_FlowCellA"],sample_expr_2 = SEQC_expr_BGI[,"A_1_L01_FlowCellA"])


(!sapply(train_and_test_res$training_res$final_model$selected.pairs.list$`3`,function(x){
  compare_two_sample(gene_order = x,
                     sample_expr_1 = SEQC_expr_AGR[,"A_1_L01_FlowCellB"],
                     sample_expr_2 = SEQC_expr_BGI[,"A_1_L01_FlowCellB"])
  
})) %>% which()


cbind(SEQC_expr_AGR[c("CRIM1","DOLPP1"),"A_1_L01_FlowCellB"],SEQC_expr_BGI[c("CRIM1","DOLPP1"),"A_1_L01_FlowCellB"])

cbind(SEQC_expr_AGR[c("TRIM3","TAF5"),"A_1_L01_FlowCellB"],SEQC_expr_BGI[c("TRIM3","TAF5"),"A_1_L01_FlowCellB"])
cbind(SEQC_expr_AGR[c("APBB2","MAML2"),"A_1_L01_FlowCellB"],SEQC_expr_BGI[c("APBB2","MAML2"),"A_1_L01_FlowCellB"])



cbind(ILM_refseq_gene_AGR[ILM_refseq_gene_AGR$Symbol%in%c("CRIM1","DOLPP1"),c("Symbol","A_1_L01_FlowCellB"),drop=FALSE],
      ILM_refseq_gene_BGI[ILM_refseq_gene_BGI$Symbol%in%c("CRIM1","DOLPP1"),c("Symbol","A_1_L01_FlowCellB"),drop=FALSE])





AIMS_diff_samples <- all_genes_predict_res[all_genes_predict_res$AIMS_predict.new!=all_genes_predict_res$AIMS_predict.old,
                                           c(1,2,8)]

sapply(AIMS_diff_samples$barcode,function(y){
  print(AIMS_diff_samples[AIMS_diff_samples$barcode==y,])
  a <- sapply(train_and_test_res$training_res$final_model$selected.pairs.list$`3`,function(x){
    compare_two_sample(gene_order = x,
                       sample_expr_1 = SEQC_expr_AGR[,y],
                       sample_expr_2 = SEQC_expr_BGI[,y])
    
  }) %>% sum
  b <- sapply(train_and_test_res$training_res$final_model$selected.pairs.list$`4`,function(x){
    compare_two_sample(gene_order = x,
                       sample_expr_1 = SEQC_expr_AGR[,y],
                       sample_expr_2 = SEQC_expr_BGI[,y])
    
  }) %>% sum
  cbind(a,b)
}) %>% t() %>% as.data.frame() %>% View





sapply(train_and_test_res$training_res$final_model$selected.pairs.list$`3`,function(x){
  expr_1_direct <- find_expr_direction(gene_order = x,sample_expr = SEQC_expr_AGR[,"A_1_L01_FlowCellB"])
  expr_1_direct == x
}) %>% sum()

num_of_consistent_direct <- function(gene_order_list,sample_expr){
  (sapply(gene_order_list,function(x){
    sample_direct <- find_expr_direction(gene_order = x,sample_expr = sample_expr)
    sample_direct == x
  })) %>% sum
}

num_of_consistent_direct(gene_order_list = train_and_test_res$training_res$final_model$selected.pairs.list$`3`,
                         sample_expr = SEQC_expr_AGR[,"A_1_L01_FlowCellB"])


sapply(AIMS_diff_samples$barcode,function(x){
  print(AIMS_diff_samples[AIMS_diff_samples$barcode==x,])
  sample_1<-num_of_consistent_direct(gene_order_list=train_and_test_res$training_res$final_model$selected.pairs.list$`3`,
                                     sample_expr = SEQC_expr_AGR[,x])
  sample_2<-num_of_consistent_direct(gene_order_list=train_and_test_res$training_res$final_model$selected.pairs.list$`3`,
                                     sample_expr = SEQC_expr_BGI[,x])
  c(sample_1,sample_2)
}) %>% t() %>% View



a <- predict.one.vs.all.tsp(D = SEQC_expr_AGR[,2,drop=FALSE],GeneName = rownames(SEQC_expr_AGR),
                            one.vs.all.tsp = train_and_test_res$training_res$final_model)

b <- predict.one.vs.all.tsp(D = SEQC_expr_BGI[,2,drop=FALSE],GeneName = rownames(SEQC_expr_BGI),
                            one.vs.all.tsp = train_and_test_res$training_res$final_model)


num_of_consistent_direct(gene_order_list = train_and_test_res$training_res$final_model$selected.pairs.list$`4`,
                         sample_expr = SEQC_expr_AGR[,2]+10)
num_of_consistent_direct(gene_order_list = train_and_test_res$training_res$final_model$selected.pairs.list$`3`,
                         sample_expr = SEQC_expr_AGR[,2]+10)


num_of_consistent_direct(gene_order_list = train_and_test_res$training_res$final_model$selected.pairs.list$`4`,
                         sample_expr = SEQC_expr_BGI[,2])






all_genes_kappa %>% as.data.frame()


all_genes_predict_res[,c(6,12)] %>% head(20)


selected_columns <- c(2,8)+5
colnames(all_genes_predict_res)[selected_columns]

all_genes_predict_res[all_genes_predict_res[,selected_columns[1]]!=all_genes_predict_res[,selected_columns[2]],
                      c(1,selected_columns)] %>% filter(grepl("^D",barcode)) %>% dim()






sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = x$PAM_res$predicted_confusions$test_set$PAM_and_real$table %>% as.matrix.data.frame)
}) %>% which.max()

ClaNC_and_PAM_train_intrinsic$`The 8st errors`$PAM_res$training_model$centroids %>% dim()


PAM_cntrd <- ClaNC_and_PAM_train_intrinsic$`The 8st errors`$PAM_res$training_model$centroids
plot(SEQC_expr_AGR[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),"A_1_L05_FlowCellB"],
     SEQC_expr_BGI[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_BGI)),"A_1_L05_FlowCellB"])

cor(SEQC_expr_AGR[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),"A_1_L05_FlowCellB"],
    SEQC_expr_BGI[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_BGI)),"A_1_L05_FlowCellB"],
    method = "spearman")


all_genes_predict_res[all_genes_predict_res$ClaNC_PAM_predict.new==all_genes_predict_res$ClaNC_PAM_predict.old,
                      c(1,6,12)] %>% head()

plot(SEQC_expr_AGR[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),"A_1_L01_FlowCellA"],
     SEQC_expr_BGI[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_BGI)),"A_1_L01_FlowCellA"])

cor(SEQC_expr_AGR[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),"A_1_L01_FlowCellA"],
    SEQC_expr_BGI[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_BGI)),"A_1_L01_FlowCellA"],
    method = "spearman")





cor(SEQC_expr_AGR[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),"A_1_L05_FlowCellB"],
    PAM_cntrd[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),],
    method = "spearman")

cor(SEQC_expr_BGI[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_BGI)),"A_1_L05_FlowCellB"],
    PAM_cntrd[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),],
    method = "spearman")


cor(SEQC_expr_AGR[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),"A_1_L01_FlowCellA"],
    PAM_cntrd[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),],
    method = "spearman")

cor(SEQC_expr_BGI[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_BGI)),"A_1_L01_FlowCellA"],
    PAM_cntrd[intersect(rownames(PAM_cntrd),rownames(SEQC_expr_AGR)),],
    method = "spearman")



#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
a <- predict_subtype(tumor_expr_dat = cbind(SEQC_expr_AGR[,"A_1_L05_FlowCellB",drop=FALSE],
                                            SEQC_expr_BGI[,"A_1_L05_FlowCellB",drop=FALSE]),
                                      training_and_testing_res = ClaNC_and_PAM_train,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))





#*******************************
#*SEQC clustering with TCGA BRCA
#*******************************
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


tumor_expr <- apply(tumor_expr,2,function(x){
  (x-mean(x))/sd(x)
})




SEQC_ILM_AGR_FPKM_entrez <- apply(ILM_refseq_gene_AGR[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_AGR[,3])

SEQC_ILM_AGR_FPKM_symbol <- apply(SEQC_ILM_AGR_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_AGR[,2])

SEQC_expr_AGR <- log2(SEQC_ILM_AGR_FPKM_symbol+1)


SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_BGI[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_BGI[,3])

SEQC_ILM_BGI_FPKM_symbol <- apply(SEQC_ILM_BGI_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_BGI[,2])

SEQC_expr_BGI <- log2(SEQC_ILM_BGI_FPKM_symbol+1)


SEQC_expr_AGR <- apply(SEQC_expr_AGR,2,function(x){
  (x-mean(x))/sd(x)
})

SEQC_expr_BGI <- apply(SEQC_expr_BGI,2,function(x){
  (x-mean(x))/sd(x)
})



dat_1 <- SEQC_expr_AGR[intersect(rownames(tumor_expr),rownames(SEQC_expr_AGR)),]
colnames(dat_1) <- paste0(colnames(dat_1),"_AGR")
dat_2 <- SEQC_expr_BGI[intersect(rownames(tumor_expr),rownames(SEQC_expr_AGR)),]
colnames(dat_2) <- paste0(colnames(dat_2),"_BGI")

TCGA_and_SEQC_data <- cbind(tumor_expr[intersect(rownames(tumor_expr),rownames(SEQC_expr_AGR)),],
                            dat_1,
                            dat_2)
sample_color <- rep(c("#E41A1C","#377EB8","#4DAF4A"),
                    times=c(ncol(tumor_expr),ncol(dat_1),ncol(dat_2)))
names(sample_color) <- colnames(TCGA_and_SEQC_data)


clustering_res <- hclust(d = dist(t(TCGA_and_SEQC_data)),method = "ward.D") %>% as.dendrogram()


pdf("../Results/Figure/Figure_5_supp_dendrogram_Z_score.pdf",width = 20,height = 10)
clustering_res %>%  
  set("labels_cex",0.1) %>% 
  #set("branches_k_color",1) %>%
  # assign_values_to_leaves_edgePar(value = ifelse(grepl("TCGA",labels(clustering_res)),"#E41A1C","#377EB8"),
  #                                 edgePar = "col") %>% 
  # set("labels_colors",value = rep(c("#E41A1C","#377EB8","#4DAF4A"),
  #                                 times=c(ncol(tumor_expr),ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI)))) %>%
  set("labels_colors",value = sample_color[labels(clustering_res)]) %>% 
  #assign_values_to_branches_edgePar(value = sample_color[labels(clustering_res)]) %>% 
  #color_branches(k=3,col = c("#E41A1C","#377EB8","#4DAF4A")) %>% 
  # hang.dendrogram(hang=5) %>% 
  plot(main="Z-Score Normalized data",cex.main=2)

legend("topright",legend = c("TCGA BRCA","SEQC AGR","SEQC BGI"),
       col=c("#E41A1C","#377EB8","#4DAF4A"),
       lwd=1,bty = "n",cex = 2)
dev.off()



SEQC_ILM_AGR_FPKM_entrez <- apply(ILM_refseq_gene_AGR[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_AGR[,3])

SEQC_ILM_AGR_FPKM_symbol <- apply(SEQC_ILM_AGR_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_AGR[,2])

SEQC_expr_AGR <- log2(SEQC_ILM_AGR_FPKM_symbol+1)


SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_BGI[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_BGI[,3])

SEQC_ILM_BGI_FPKM_symbol <- apply(SEQC_ILM_BGI_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_BGI[,2])

SEQC_expr_BGI <- log2(SEQC_ILM_BGI_FPKM_symbol+1)


SEQC_ILM_CNL_FPKM_entrez <- apply(ILM_refseq_gene_CNL[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_CNL[,3])

SEQC_ILM_CNL_FPKM_symbol <- apply(SEQC_ILM_CNL_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_CNL[,2])

SEQC_expr_CNL <- log2(SEQC_ILM_CNL_FPKM_symbol+1)





dat_1 <- SEQC_expr_AGR[rownames(SEQC_expr_AGR),]
colnames(dat_1) <- paste0(colnames(dat_1),"_AGR")
dat_2 <- SEQC_expr_BGI[rownames(SEQC_expr_AGR),]
colnames(dat_2) <- paste0(colnames(dat_2),"_BGI")
dat_2 <- dat_2[,!grepl("_5_",colnames(dat_2))]


TCGA_and_SEQC_data <- cbind(dat_1,dat_2)
sample_color <- rep(c("#377EB8","#4DAF4A"),
                    times=c(ncol(dat_1),ncol(dat_2)))
names(sample_color) <- colnames(TCGA_and_SEQC_data)

clustering_res <- hclust(d = dist(t(TCGA_and_SEQC_data)),method = "ward.D") %>% as.dendrogram()

pdf("../Results/Figure/Figure_5_supp_dendrogram_SEQC_original_no_sample_5.pdf",width = 20,height = 10)
clustering_res %>%  
  set("labels_cex",0.1) %>% 
  #set("branches_k_color",1) %>%
  # assign_values_to_leaves_edgePar(value = ifelse(grepl("TCGA",labels(clustering_res)),"#E41A1C","#377EB8"),
  #                                 edgePar = "col") %>% 
  # set("labels_colors",value = rep(c("#E41A1C","#377EB8","#4DAF4A"),
  #                                 times=c(ncol(tumor_expr),ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI)))) %>%
  set("labels_colors",value = sample_color[labels(clustering_res)]) %>% 
  #assign_values_to_branches_edgePar(value = sample_color[labels(clustering_res)]) %>% 
  #color_branches(k=3,col = c("#E41A1C","#377EB8","#4DAF4A")) %>% 
  # hang.dendrogram(hang=5) %>% 
  plot(main="Z-Score Normalized data",cex.main=2)

legend("topright",legend = c("SEQC AGR","SEQC BGI"),
       col=c("#377EB8","#4DAF4A"),
       lwd=1,bty = "n",cex = 2)
dev.off()



TCGA_and_SEQC_data <- cbind(dat_1)
color_lab <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
names(color_lab) <- c("A","B","C","D")
sample_color <- color_lab[colnames(dat_1) %>% substr(1,1)]
names(sample_color) <- colnames(TCGA_and_SEQC_data)

clustering_res <- hclust(d = dist(t(TCGA_and_SEQC_data)),method = "ward.D") %>% as.dendrogram()

pdf("../Results/Figure/Figure_5_supp_dendrogram_SEQC_AGR_original.pdf",width = 20,height = 10)
clustering_res %>%  
  set("labels_cex",0.5) %>% 
  #set("branches_k_color",1) %>%
  # assign_values_to_leaves_edgePar(value = ifelse(grepl("TCGA",labels(clustering_res)),"#E41A1C","#377EB8"),
  #                                 edgePar = "col") %>% 
  # set("labels_colors",value = rep(c("#E41A1C","#377EB8","#4DAF4A"),
  #                                 times=c(ncol(tumor_expr),ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI)))) %>%
  set("labels_colors",value = sample_color[labels(clustering_res)]) %>% 
  #assign_values_to_branches_edgePar(value = sample_color[labels(clustering_res)]) %>% 
  #color_branches(k=3,col = c("#E41A1C","#377EB8","#4DAF4A")) %>% 
  # hang.dendrogram(hang=5) %>% 
  plot(main="Z-Score Normalized data",cex.main=2)

legend("topright",legend = c("SEQC AGR","SEQC BGI"),
       col=c("#377EB8","#4DAF4A"),
       lwd=1,bty = "n",cex = 2)
dev.off()




dat_1 <- SEQC_expr_CNL[rownames(SEQC_expr_CNL),]
colnames(dat_1) <- paste0(colnames(dat_1),"_CNL")
dat_2 <- SEQC_expr_BGI[rownames(SEQC_expr_BGI),]
colnames(dat_2) <- paste0(colnames(dat_2),"_BGI")

TCGA_and_SEQC_data <- cbind(dat_1,dat_2)
sample_color <- rep(c("#377EB8","#4DAF4A"),
                    times=c(ncol(dat_1),ncol(dat_2)))
names(sample_color) <- colnames(TCGA_and_SEQC_data)

clustering_res <- hclust(d = dist(t(TCGA_and_SEQC_data)),method = "ward.D") %>% as.dendrogram()

pdf("../Results/Figure/Figure_5_supp_dendrogram_SEQC_BGI_CNL_original_2.pdf",width = 20,height = 10)
clustering_res %>%  
  set("labels_cex",0.1) %>% 
  #set("branches_k_color",1) %>%
  assign_values_to_leaves_edgePar(value = ifelse(grepl("TCGA",labels(clustering_res)),"#E41A1C","#377EB8"),
                                  edgePar = "col") %>%
  # set("labels_colors",value = rep(c("#E41A1C","#377EB8","#4DAF4A"),
  #                                 times=c(ncol(tumor_expr),ncol(SEQC_expr_AGR),ncol(SEQC_expr_BGI)))) %>%
  set("labels_colors",value = sample_color[labels(clustering_res)]) %>% 
  #assign_values_to_branches_edgePar(value = sample_color[labels(clustering_res)]) %>% 
  #color_branches(k=3,col = c("#E41A1C","#377EB8","#4DAF4A")) %>% 
  # hang.dendrogram(hang=5) %>% 
  # plot(main="Z-Score Normalized data",cex.main=2)
  plot(main="Original Expression data",cex.main=2)

legend("topright",legend = c("SEQC BGI","SEQC CNL"),
       col=c("#377EB8","#4DAF4A"),
       lwd=1,bty = "n",cex = 2)
dev.off()



expr_dt <- TCGA_and_SEQC_data
pc.cr<-prcomp(t(expr_dt),retx = TRUE)
pc<-round(summary(pc.cr)$importance[2,],2)

palette(c(RColorBrewer::brewer.pal(9,"Set1")))
plot(pc.cr$x[,1:2],xlab=paste("PC1",pc[1]),
     ylab=paste("PC2",pc[2]),col=rep(c(1,2),times=c(ncol(dat_1),ncol(dat_2))),
     pch=rownames(pc.cr$x) %>% substr(1,1) %>% as.factor() %>% as.integer())
legend("topleft",legend = unique(str_remove_all(colnames(expr_dat),".*_")),col=1:2,
       pch=16,bty="n")

legend("topright",legend = unique(rownames(pc.cr$x) %>% substr(1,1)),
       pch=rownames(pc.cr$x) %>% substr(1,1) %>% unique() %>% as.factor() %>% as.integer(),bty="n")




