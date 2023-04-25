##classification consistency between RNA-seq and microarray
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sigclust2)
library(tidyverse)
library(dendextend)
library(pheatmap)
library(pamr)
library(survminer)
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
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))

#AIMS intersected intrinsic genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))

#all intersect genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
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

#instersected intrinsic genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData")
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


#RNA-seq predicting results for test sample of each predicting strategy







#*******************************
#*array data model predict array
#*******************************
#AIMS all intersect genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_AIMS_training_res.RData")
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))

#AIMS intersected intrinsic genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))



#all intersect genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
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



#instersected intrinsic genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData")
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




#*****************************************
#*array and RNA-seq predict same sample***
#*using array and RNA-seq data repectively
#*****************************************
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Results/")

#AIMS all intersect genes
RNA_seq_intersect_AIMS<-obtain_loaded_dat("TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")
array_intersect_AIMS<-obtain_loaded_dat("TCGA_array_intersect_gene_consensus_sample_AIMS_training_res.RData")
kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_AIMS,array_pred = array_intersect_AIMS)$overall["Kappa"]

#AIMS intrinsic genes
RNA_seq_intr<-obtain_loaded_dat("TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
array_intr <- obtain_loaded_dat("TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intr,array_pred = array_intr)$overall["Kappa"]

##all intersect genes
RNA_seq_intersect <- obtain_loaded_dat("TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
array_intersect <- obtain_loaded_dat("TCGA_array_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
#PAM predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})
#PAM + Cor predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})

#ClaNC predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$ClaNC_res$testing_predicted_res")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})
#ClaNC + PAM predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})
#ClaNC + PAM + Cor predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})


##all intrinsic genes

RNA_seq_intersect_intr <- obtain_loaded_dat(paste0("TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene",
                                                   "_ClaNC_and_PAM_training_res.RData"))
array_intersect_intr <- obtain_loaded_dat(paste0("TCGA_array_intersect_gene_consensus_sample_intrinsic_gene",
                                                 "_ClaNC_and_PAM_training_res.RData"))

#PAM predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
})
#PAM + Cor predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
}) %>% mean()

#ClaNC predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$ClaNC_res$testing_predicted_res")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
}) %>% mean()
#ClaNC + PAM predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
}) %>% mean()
#ClaNC + PAM + Cor predicting kappa
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
}) %>% mean()





# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

load("/mnt/Miscrosoft/TCGA_project/Data/R_data/TCGA-BRCA/Expr_and_pheno.RData")
expr_dat <- X1_genename
expr_dat <- log2(expr_dat+1)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]

tumor_sample_ID_dat <- A %>% 
  dplyr::filter(sample_type=="Primary Tumor") %>% dplyr::select(barcode,sample_type)


#tumor sample expr
tumor_expr <- expr_dat[,tumor_sample_ID_dat$barcode]


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
  pamr_centroids <- pamr_centroids[pamr_survived_genes,]
  pamr_centroids_test_cor <- cor(x = expr_dat[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
  test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
}


#predict subtypes of each sample
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
           PAM_predict <- pamr.predict(fit = PAM_model$training_model,
                                       newx = tumor_expr_dat[rownames(PAM_model$training_model$centroids),],
                                       threshold = PAM_model$delt)
           names(PAM_predict) <- colnames(tumor_expr_dat)
           
           #Select PAM + Cor model with highest MCC
           predict_highest_MCC_Cor <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
           }))
           PAM_model <- training_and_testing_res[[predict_highest_MCC_Cor]]
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
           PAM_predict <- pamr.predict(fit = PAM_model$training_model,
                                       newx = tumor_expr_dat[rownames(PAM_model$training_model$centroids),],
                                       threshold = PAM_model$delt)
           names(PAM_predict) <- colnames(tumor_expr_dat)
           
           #Select PAM + Cor model with highest MCC
           predict_highest_MCC_Cor <- which.max(sapply(training_and_testing_res,function(x){
             mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
           }))
           PAM_model <- training_and_testing_res[[predict_highest_MCC_Cor]]$PAM_res
           PAM_Cor_predict <- cor_predict(expr_dat = tumor_expr_dat,PAM_model = PAM_model)
           
           cbind(ClaNC_predict,ClaNC_PAM_predict=PAM_predict,ClaNC_PAM_Cor_predict=PAM_Cor_predict) %>% 
             as.data.frame %>% rownames_to_column("barcode")
         })
  
}



# tumor_expr_dat <- tumor_expr
# training_and_testing_res <- train_and_test_res_intrinsic
# training_and_testing_res <- PAM_and_PAM_plus_Cor
# training_and_testing_res <- ClaNC_and_PAM_train
# load("../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")
# class_label <- levels(as.factor(datset_for_model$y))


#RNA-seq All intersect genes predicting res
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
load("../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")
load("TCGA_BRCA_clinical_dat.RData")

#AIMS
all_gene_AIMS_res <- predict_subtype(tumor_expr_dat = tumor_expr,
                                     training_and_testing_res = train_and_test_res,
                                     predicting_type = "AIMS")

#PAM and PAM + Cor
all_gene_PAM_res <- predict_subtype(tumor_expr_dat = tumor_expr,training_and_testing_res = PAM_and_PAM_plus_Cor,
                                    predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = tumor_expr,training_and_testing_res = ClaNC_and_PAM_train,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))

#combined_res
combined_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res),left_join,by="barcode")


dat_of_clinical <- combined_res %>% left_join(A,by="barcode") %>% 
  select(1:7,barcode,patient,
         days_to_last_follow_up) %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>%
  select(1:7,barcode,patient,
         days_to_last_follow_up,death_days_to,vital_status) %>%
  apply(2,as.character) %>% as.data.frame() %>%
  mutate(days_to_death=str_replace_all(death_days_to,'[^0-9]','')) %>%
  mutate(days_to_last_follow_up=ifelse(days_to_death=='',days_to_last_follow_up,days_to_death)) %>%
  mutate(days_to_death=ifelse(days_to_death=='',NA,as.integer(days_to_death))) %>%
  mutate(days_to_last_follow_up=as.integer(days_to_last_follow_up)) %>% filter(!is.na(vital_status))

surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))

fit <- survfit(Surv(days_to_last_follow_up,vital_status)~AIMS_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)

# sapply(colnames(surv_dat)[2:7],function(x){
#   fit <- survfit(as.formula(paste0('Surv(days_to_last_follow_up,vital_status)~',eval(x))),data = surv_dat)
#   fit <- survfit(reformulate(x,'Surv(days_to_last_follow_up,vital_status)'),data = surv_dat)
#   coxph(reformulate(x,'Surv(days_to_last_follow_up,vital_status)'),data = surv_dat)
#   
#   #g_surv <- ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
#   #                     font.tickslab = 15,termlabels = as.character(1:5))
# })
# x <- "AIMS_predict"
# fit <- survfit(Surv(surv_dat$days_to_last_follow_up,surv_dat$vital_status)~surv_dat[[x]],data = surv_dat)
  
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~PAM_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~PAM_Cor_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~ClaNC_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~ClaNC_PAM_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~ClaNC_PAM_Cor_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)









#RNA-seq Intrinsic genes predicting res

#RNA-seq All intersect genes predicting res
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData")
load("../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")
load("TCGA_BRCA_clinical_dat.RData")

#AIMS
all_gene_AIMS_res <- predict_subtype(tumor_expr_dat = tumor_expr,
                                     training_and_testing_res = train_and_test_res_intrinsic,
                                     predicting_type = "AIMS")

#PAM and PAM + Cor
all_gene_PAM_res <- predict_subtype(tumor_expr_dat = tumor_expr,
                                    training_and_testing_res = PAM_and_PAM_plus_Cor_intrinsic,
                                    predicting_type = "PAM")

#ClaNC, ClaNC + PAM and ClaNC + PAM + Cor
all_gene_ClaNC_res <- predict_subtype(tumor_expr_dat = tumor_expr,
                                      training_and_testing_res = ClaNC_and_PAM_train_intrinsic,
                                      predicting_type = "ClaNC",class_label = levels(as.factor(datset_for_model$y)))

#combined_res
combined_res <- reduce(list(all_gene_AIMS_res,all_gene_PAM_res,all_gene_ClaNC_res),left_join,by="barcode")


dat_of_clinical <- combined_res %>% left_join(A,by="barcode") %>% 
  select(1:7,barcode,patient,
         days_to_last_follow_up) %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>%
  select(1:7,barcode,patient,
         days_to_last_follow_up,death_days_to,vital_status) %>%
  apply(2,as.character) %>% as.data.frame() %>%
  mutate(days_to_death=str_replace_all(death_days_to,'[^0-9]','')) %>%
  mutate(days_to_last_follow_up=ifelse(days_to_death=='',days_to_last_follow_up,days_to_death)) %>%
  mutate(days_to_death=ifelse(days_to_death=='',NA,as.integer(days_to_death))) %>%
  mutate(days_to_last_follow_up=as.integer(days_to_last_follow_up)) %>% filter(!is.na(vital_status))

surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))

fit <- survfit(Surv(days_to_last_follow_up,vital_status)~AIMS_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)


fit <- survfit(Surv(days_to_last_follow_up,vital_status)~PAM_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~PAM_Cor_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~ClaNC_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~ClaNC_PAM_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~ClaNC_PAM_Cor_predict,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
           font.tickslab = 15)











