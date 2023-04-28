#****************************************************************************
#*Table 2 and Supplementary Table 2 in manuscript
#****************************************************************************

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



setwd("./Scripts/")
source("./clanc.R")
source("./all_functions_20220302.R")
source("./trainAIMS_2.R")


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
setwd("../Data/")



#**************************
#*Supplementary Table 2****
#**************************

#AIMS all intersect genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")

#****************************
#*RNA-seq Gene-pairs all gene
#****************************
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))

#AIMS intersected intrinsic genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")

#**********************************
#*RNA-seq Gene-pairs intrinsic gene
#**********************************
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))

#all intersect genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")

#*********************
#*RNA-seq PAM all gene
#*********************
mean(sapply(PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))

#******************************
#*RNA-seq PAM+Spearman all gene
#******************************
mean(sapply(PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))

#***********************
#*RNA-seq ClaNC all gene
#***********************
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))

#***************************
#*RNA-seq ClaNC+PAM all gene
#***************************
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))

#************************************
#*RNA-seq ClaNC+PAM+Spearman all gene
#************************************
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))

#instersected intrinsic genes
load("../Results/TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData")

#***************************
#*RNA-seq PAM intrinsic gene
#***************************
mean(sapply(PAM_and_PAM_plus_Cor_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))

#************************************
#*RNA-seq PAM+Spearman intrinsic gene
#************************************
mean(sapply(PAM_and_PAM_plus_Cor_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))

#*****************************
#*RNA-seq ClaNC intrinsic gene
#*****************************
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))

#*********************************
#*RNA-seq ClaNC+PAM intrinsic gene
#*********************************
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))

#******************************************
#*RNA-seq ClaNC+PAM+Spearman intrinsic gene
#******************************************
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))


#RNA-seq predicting results for test sample of each predicting strategy







#*******************************
#*array data model predict array
#*******************************
#AIMS all intersect genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_AIMS_training_res.RData")

#**************************
#*Array Gene-pairs all gene
#**************************
train_and_test_res$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res$test_confusion))

#AIMS intersected intrinsic genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
#********************************
#*Array Gene-pairs intrinsic gene
#********************************
train_and_test_res_intrinsic$test_confusion
mltools::mcc(confusionM = as.matrix.data.frame(train_and_test_res_intrinsic$test_confusion))



#all intersect genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")

#*******************
#*Array PAM all gene
#*******************
mean(sapply(PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))

#****************************
#*Array PAM+Spearman all gene
#****************************
mean(sapply(PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))

#*********************
#*Array ClaNC all gene
#*********************
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))

#*************************
#*Array ClaNC+PAM all gene
#*************************
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))

#**********************************
#*Array ClaNC+PAM+Spearman all gene
#**********************************
mean(sapply(ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))



#instersected intrinsic genes
load("../Results/TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_ClaNC_and_PAM_training_res.RData")

#*************************
#*Array PAM intrinsic gene
#*************************
mean(sapply(PAM_and_PAM_plus_Cor_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))

#**********************************
#*Array PAM+Spearman intrinsic gene
#**********************************
mean(sapply(PAM_and_PAM_plus_Cor_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))

#***************************
#*Array ClaNC intrinsic gene
#***************************
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))

#*******************************
#*Array ClaNC+PAM intrinsic gene
#*******************************
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))

#****************************************
#*Array ClaNC+PAM+Spearman intrinsic gene
#****************************************
mean(sapply(ClaNC_and_PAM_train_intrinsic,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))




#*****************************************
#*array and RNA-seq predict same sample***
#*using array and RNA-seq data repectively
#*****************************************
setwd("../Results/")


#************
#*Table 2****
#************

#AIMS all intersect genes
RNA_seq_intersect_AIMS<-obtain_loaded_dat("TCGA_RNA_seq_intersect_gene_consensus_sample_AIMS_training_res.RData")
array_intersect_AIMS<-obtain_loaded_dat("TCGA_array_intersect_gene_consensus_sample_AIMS_training_res.RData")

#********************
#*Gene-pairs all gene
#********************
kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_AIMS,array_pred = array_intersect_AIMS)$overall["Kappa"]

#AIMS intrinsic genes
RNA_seq_intr<-obtain_loaded_dat("TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")
array_intr <- obtain_loaded_dat("TCGA_array_intersect_gene_consensus_sample_intrinsic_gene_AIMS_training_res.RData")

#**************************
#*Gene-pairs intrinsic gene
#**************************
kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intr,array_pred = array_intr)$overall["Kappa"]



##all intersect genes
RNA_seq_intersect <- obtain_loaded_dat("TCGA_RNA_seq_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")
array_intersect <- obtain_loaded_dat("TCGA_array_intersect_gene_consensus_sample_ClaNC_and_PAM_training_res.RData")

#*************
#*PAM all gene
#*************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect$PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))]

#**********************
#*PAM+Spearman all gene
#**********************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect$PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))]

#***************
#*ClaNC all gene
#***************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$ClaNC_res$testing_predicted_res")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect$ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))]

#*******************
#*ClaNC+PAM all gene
#*******************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect$ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))]

#****************************
#*ClaNC+PAM+Spearman all gene
#****************************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect,array_pred = array_intersect,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect$ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))]




##all intrinsic genes
RNA_seq_intersect_intr <- obtain_loaded_dat(paste0("TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene",
                                                   "_ClaNC_and_PAM_training_res.RData"))
array_intersect_intr <- obtain_loaded_dat(paste0("TCGA_array_intersect_gene_consensus_sample_intrinsic_gene",
                                                 "_ClaNC_and_PAM_training_res.RData"))


#*******************
#*PAM intrinsic gene
#*******************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect_intr$PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$PAM_and_real$table))
}))]

#****************************
#*PAM+Spearman intrinsic gene
#****************************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$PAM_and_PAM_plus_Cor$`",nth,"`$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect_intr$PAM_and_PAM_plus_Cor,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$predicted_confusions$test_set$Cor_and_real$table))
}))]

#*********************
#*ClaNC intrinsic gene
#*********************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$ClaNC_res$testing_predicted_res")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect_intr$ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$ClaNC_res$predicted_confusions$test_set$table))
}))]

#*************************
#*ClaNC+PAM intrinsic gene
#*************************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$pam_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect_intr$ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$PAM_and_real$table))
}))]

#**********************************
#*ClaNC+PAM+Spearman intrinsic gene
#**********************************
sapply(1:20,function(x){
  nth <- paste0("The ",x,"st errors")
  compare_part <- paste0("$ClaNC_and_PAM_train$`",nth,"`$PAM_res$predict_res$test_set$cor_predict")
  kappa_between_array_and_RNA_seq(RNA_pred = RNA_seq_intersect_intr,array_pred = array_intersect_intr,
                                  compare_part=compare_part)$overall["Kappa"]
})[which.max(sapply(RNA_seq_intersect_intr$ClaNC_and_PAM_train,function(x){
  mltools::mcc(confusionM = as.matrix.data.frame(x$PAM_res$predicted_confusions$test_set$Cor_and_real$table))
}))]















