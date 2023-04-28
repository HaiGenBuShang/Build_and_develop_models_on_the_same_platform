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
RNA_seq_expr_for_clust <- tumor_intrinsic_expr
RNA_seq_expr <- tumor_expr


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
tumor_expr <- tumor_expr[!apply(tumor_expr,1,function(x)any(is.na(x))),]


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
array_expr <- tumor_expr


# For clustering intrinsic gene clustering --------------------------------

RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[,
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,scale))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,scale))

combined_dat <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_RNAseq"),
                            paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_array"))


#clustering methods were changed from ward.D to Ward.D
clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D")
gene_clustering_res <- hclust(d = dist(combined_dat),method = "ward.D")



samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1) %>% mutate(cluster_num=as.character(cluster_num))


consensus_samples_array_and_RNA_seq <- samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>%
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

# samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
#   group_by(samples,cluster_num) %>% count() %>% filter(n!=2)





samples_of_5_classes <- cutree(clustering_res,k=5) %>% as.matrix() %>% as.data.frame() %>% 
  rownames_to_column("barcode") %>% rename(cluster_num=2) %>% mutate(cluster_num=as.character(cluster_num)) %>% 
  mutate(samples=str_replace_all(barcode,"_.*",""))

consensus_samples_of_5_classes <- samples_of_5_classes %>% 
  inner_join(consensus_samples_array_and_RNA_seq,by="samples") %>% 
  mutate(data_type=str_replace_all(barcode,".*_","")) %>% 
  pivot_wider(id_cols = -barcode,names_from = data_type,values_from = data_type)



# train model using AIMS idea ---------------------------------------------

#intersect genes
intersect_intrinsic_genes <- intersect(intersect_genes,intrinsic_genes[,1])
intersect_all_genes <- intersect(rownames(RNA_seq_expr),rownames(array_expr))

#*************
#*RNA-seq data
#*************
tumor_expr <- RNA_seq_expr[intersect_all_genes,consensus_samples_of_5_classes$samples]
cluster_sample_subtypes <- consensus_samples_of_5_classes$cluster_num.x
tumor_intrinsic_expr <- RNA_seq_expr[intersect_intrinsic_genes,consensus_samples_of_5_classes$samples]


dat_for_model <- tumor_expr
datset_for_model <- list(x=dat_for_model,y=cluster_sample_subtypes)
datset_for_model_intrinsic_g <- list(x=tumor_intrinsic_expr,y=cluster_sample_subtypes)
save(datset_for_model,datset_for_model_intrinsic_g,
    file = "../Results/datset_for_model_and_dataset_for_model_intrinsic_g.RData")

res_prefix <- "TCGA_RNA_seq_intersect_gene_consensus_sample"
res_prefix_intrinsic <- "TCGA_RNA_seq_intersect_gene_consensus_sample_intrinsic_gene"





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

# #**************
# #*Script Part 1
# #**************
# set.seed(12345678)
# train_and_test_with_label <- produce_train_test_set(expr_with_label = datset_for_model)
# train_and_test_res <- AIMS_train_and_pred(train_and_test_data_with_lable = train_and_test_with_label,
#                                           PREFIX = res_prefix,
#                                           k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res,file = paste0(res_prefix,"_AIMS_training_res.RData"))
# message("Finished all genes AIMS!")

# #**************
# #*Script Part 2
# #**************
# set.seed(12345678)
# #using intrinsic genes
# train_and_test_with_label_intrinsic <- produce_train_test_set(expr_with_label = datset_for_model_intrinsic_g)
# train_and_test_res_intrinsic<-AIMS_train_and_pred(
#   train_and_test_data_with_lable=train_and_test_with_label_intrinsic,PREFIX = res_prefix_intrinsic,
#   k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res_intrinsic,file = paste0(res_prefix_intrinsic,"_AIMS_training_res.RData"))
# message("Finished intrinsic genes AIMS!")


# #**************
# #*Script Part 3
# #**************
# set.seed(12345678)
# # train model using PAM50 strategies
# # using all genes
# PAM_and_PAM_plus_Cor <- n_times_compare(expr_with_label = datset_for_model,n_times = 20,train_proportion = 2/3,
#                                         auto_pam_delt = TRUE,n_threshold = 50,
#                                         auto_delt_by="min_CV_error",prior="class")
# 
# ClaNC_and_PAM_train <- n_times_CLaNC_and_PAM_com(expr_with_label = datset_for_model,ntimes = 20,
#                                                  train_proportion = 2/3,
#                                                  already_train_and_test_set=FALSE,show_message = FALSE,
#                                                  prior = "class",CV_gene_number = 1:50,auto_active_genes = TRUE,
#                                                  auto_pam_delt=FALSE,manual_delt=0)
# save(PAM_and_PAM_plus_Cor,ClaNC_and_PAM_train,file = paste0(res_prefix,"_ClaNC_and_PAM_training_res.RData"))
# 
# 
# #**************
# #*Script Part 4
# #**************
# set.seed(12345678)
# #train model using PAM50 strategies
# #using intrinsic gens
# PAM_and_PAM_plus_Cor_intrinsic <- n_times_compare(expr_with_label = datset_for_model_intrinsic_g,
#                                                   n_times = 20,train_proportion = 2/3,
#                                                   auto_pam_delt = TRUE,n_threshold = 50,
#                                                   auto_delt_by="min_CV_error",prior="class")
# 
# ClaNC_and_PAM_train_intrinsic <- n_times_CLaNC_and_PAM_com(expr_with_label = datset_for_model_intrinsic_g,
#                                                            ntimes = 20,
#                                                            train_proportion = 2/3,
#                                                            already_train_and_test_set=FALSE,show_message = FALSE,
#                                                            prior = "class",CV_gene_number = 1:50,
#                                                            auto_active_genes = TRUE,
#                                                            auto_pam_delt=FALSE,manual_delt=0)
# save(PAM_and_PAM_plus_Cor_intrinsic,ClaNC_and_PAM_train_intrinsic,
#      file = paste0(res_prefix_intrinsic,"_ClaNC_and_PAM_training_res.RData"))
# 


#***********
#*array data
#***********
tumor_expr <- array_expr[intersect_all_genes,consensus_samples_of_5_classes$samples]
cluster_sample_subtypes <- consensus_samples_of_5_classes$cluster_num.x
tumor_intrinsic_expr <- array_expr[intersect_intrinsic_genes,consensus_samples_of_5_classes$samples]


dat_for_model <- tumor_expr
datset_for_model <- list(x=dat_for_model,y=cluster_sample_subtypes)
datset_for_model_intrinsic_g <- list(x=tumor_intrinsic_expr,y=cluster_sample_subtypes)
res_prefix <- "TCGA_array_intersect_gene_consensus_sample"
res_prefix_intrinsic <- "TCGA_array_intersect_gene_consensus_sample_intrinsic_gene"




# #**************
# #*Script Part 5
# #**************
# set.seed(12345678)
# train_and_test_with_label <- produce_train_test_set(expr_with_label = datset_for_model)
# train_and_test_res <- AIMS_train_and_pred(train_and_test_data_with_lable = train_and_test_with_label,
#                                           PREFIX = res_prefix,
#                                           k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res,file = paste0(res_prefix,"_AIMS_training_res.RData"))
# message("Finished all genes AIMS!")

# #**************
# #*Script Part 6
# #**************
# set.seed(12345678)
# #using intrinsic genes
# train_and_test_with_label_intrinsic <- produce_train_test_set(expr_with_label = datset_for_model_intrinsic_g)
# train_and_test_res_intrinsic<-AIMS_train_and_pred(
#   train_and_test_data_with_lable=train_and_test_with_label_intrinsic,PREFIX = res_prefix_intrinsic,
#   k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res_intrinsic,file = paste0(res_prefix_intrinsic,"_AIMS_training_res.RData"))
# message("Finished intrinsic genes AIMS!")


#**************
#*Script Part 7
#**************
set.seed(12345678)
# train model using PAM50 strategies
# using all genes
PAM_and_PAM_plus_Cor <- n_times_compare(expr_with_label = datset_for_model,n_times = 20,train_proportion = 2/3,
                                        auto_pam_delt = TRUE,n_threshold = 50,
                                        auto_delt_by="min_CV_error",prior="class")

ClaNC_and_PAM_train <- n_times_CLaNC_and_PAM_com(expr_with_label = datset_for_model,ntimes = 20,
                                                 train_proportion = 2/3,
                                                 already_train_and_test_set=FALSE,show_message = FALSE,
                                                 prior = "class",CV_gene_number = 1:50,auto_active_genes = TRUE,
                                                 auto_pam_delt=FALSE,manual_delt=0)
save(PAM_and_PAM_plus_Cor,ClaNC_and_PAM_train,file = paste0(res_prefix,"_ClaNC_and_PAM_training_res.RData"))


#**************
#*Script Part 8
#**************
set.seed(12345678)
#train model using PAM50 strategies
#using intrinsic gens
PAM_and_PAM_plus_Cor_intrinsic <- n_times_compare(expr_with_label = datset_for_model_intrinsic_g,
                                                  n_times = 20,train_proportion = 2/3,
                                                  auto_pam_delt = TRUE,n_threshold = 50,
                                                  auto_delt_by="min_CV_error",prior="class")

ClaNC_and_PAM_train_intrinsic <- n_times_CLaNC_and_PAM_com(expr_with_label = datset_for_model_intrinsic_g,
                                                           ntimes = 20,
                                                           train_proportion = 2/3,
                                                           already_train_and_test_set=FALSE,show_message = FALSE,
                                                           prior = "class",CV_gene_number = 1:50,
                                                           auto_active_genes = TRUE,
                                                           auto_pam_delt=FALSE,manual_delt=0)
save(PAM_and_PAM_plus_Cor_intrinsic,ClaNC_and_PAM_train_intrinsic,
     file = paste0(res_prefix_intrinsic,"_ClaNC_and_PAM_training_res.RData"))






