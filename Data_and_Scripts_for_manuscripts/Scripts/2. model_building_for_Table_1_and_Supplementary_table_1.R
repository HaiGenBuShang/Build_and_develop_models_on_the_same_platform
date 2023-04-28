library(TCGAbiolinks)
library(SummarizedExperiment)
library(sigclust2)
library(tidyverse)
library(dendextend)
library(pheatmap)
library(pamr)


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






set.seed(12345678)
#down load TCGA BRCA expressions
setwd("../Data/")

# array data
load("20221216_BRCA_microarray_expression.RData")

expr_dat <- BRCA_array_expr %>% apply(2,as.numeric)
rownames(expr_dat) <- rownames(BRCA_array_expr)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]

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
expr_for_clust <- tumor_intrinsic_expr

#intrinsic_genes_for_corresponding
#would be changed according to PNAS 2001, 2003,
#and myabe other genes in other papers should be integrated.
intrinsic_genes_for_corresponding <- data.frame(intrinsic_genes=rownames(expr_for_clust),
                                                stringsAsFactors = FALSE)
intrinsic_gene_dat <- intrinsic_genes_for_corresponding %>% 
  mutate(determined_subtype=sample(c("Basal","LumA","LumB","Her2","Normal"),
                                   replace = TRUE,size = nrow(intrinsic_genes_for_corresponding)))

#do clust
#use euclidean and ward.D
clust_res <- hclust(d = dist(t(expr_for_clust),method = "euclidean"),method = "ward.D")
samples_of_classes <- cutree(clust_res,k=5)

#cut to 5 classes
samples_of_classes <- cutree(clust_res,k = 5)
#plot(color_branches(clust_res,k=5),leaflab="none")

gene_clust <- hclust(d = dist(expr_for_clust,method = "euclidean"),method = "ward.D")


cluster_sample_subtypes <- as.character(samples_of_classes)


#genes intersected with RNA_seq data
load("./R_data/BRCA/Expr_and_pheno.RData")
expr_dat <- X1_genename
RNA_seq_genes <- rownames(expr_dat)
intersect_genes <- intersect(RNA_seq_genes,rownames(tumor_expr))
intersect_intrinsic_genes <- intersect(intersect_genes,intrinsic_genes[,1])


#train model using AIMS idea
dat_for_model <- tumor_expr[intersect_genes,]
datset_for_model <- list(x=dat_for_model,y=cluster_sample_subtypes)
datset_for_model_intrinsic_g <- list(x=tumor_intrinsic_expr[intersect_intrinsic_genes,],y=cluster_sample_subtypes)
res_prefix <- "TCGA_array_intersect_gene_train_and_test"
res_prefix_intrinsic <- "TCGA_array_intersect_gene_train_and_test_intrinsic_gene"

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



#**************
#*Script Part 2
#**************
set.seed(12345678)
#using intrinsic genes
train_and_test_with_label_intrinsic <- produce_train_test_set(expr_with_label = datset_for_model_intrinsic_g)
train_and_test_res_intrinsic<-AIMS_train_and_pred(
  train_and_test_data_with_lable=train_and_test_with_label_intrinsic,PREFIX = res_prefix_intrinsic,
  k.fold = 10,num.of.rules = seq(1,50,1))
save(train_and_test_res_intrinsic,file = paste0(res_prefix_intrinsic,"_AIMS_training_res.RData"))
message("Finished intrinsic genes AIMS!")



#**************
#*Script Part 3
#**************
set.seed(12345678)
# using all genes
PAM_and_PAM_plus_Cor <- n_times_compare(expr_with_label = datset_for_model,n_times = 20,train_proportion = 2/3,
                                        auto_pam_delt = TRUE,n_threshold = 50,
                                        auto_delt_by="min_CV_error",prior="class")
set.seed(12345678)
ClaNC_and_PAM_train <- n_times_CLaNC_and_PAM_com(expr_with_label = datset_for_model,ntimes = 20,
                                                 train_proportion = 2/3,
                                                 already_train_and_test_set=FALSE,show_message = FALSE,
                                                 prior = "class",CV_gene_number = 1:50,auto_active_genes = TRUE,
                                                 auto_pam_delt=FALSE,manual_delt=0)
save(PAM_and_PAM_plus_Cor,ClaNC_and_PAM_train,file = paste0(res_prefix,"_ClaNC_and_PAM_training_res.RData"))



#**************
#*Script Part 4
#**************
set.seed(12345678)
#using intrinsic gens
PAM_and_PAM_plus_Cor_intrinsic <- n_times_compare(expr_with_label = datset_for_model_intrinsic_g,
                                                  n_times = 20,train_proportion = 2/3,
                                                  auto_pam_delt = TRUE,n_threshold = 50,
                                                  auto_delt_by="min_CV_error",prior="class")
set.seed(12345678)
ClaNC_and_PAM_train_intrinsic <- n_times_CLaNC_and_PAM_com(expr_with_label = datset_for_model_intrinsic_g,
                                                           ntimes = 20,
                                                           train_proportion = 2/3,
                                                           already_train_and_test_set=FALSE,show_message = FALSE,
                                                           prior = "class",CV_gene_number = 1:50,
                                                           auto_active_genes = TRUE,
                                                           auto_pam_delt=FALSE,manual_delt=0)
save(PAM_and_PAM_plus_Cor_intrinsic,ClaNC_and_PAM_train_intrinsic,
     file = paste0(res_prefix_intrinsic,"_ClaNC_and_PAM_training_res.RData"))
