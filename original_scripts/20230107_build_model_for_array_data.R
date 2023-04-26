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


setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")


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
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

# ##try TCGA-BRCA lagacy data 
# #normalized results
# query <- GDCquery(
#   project = "TCGA-BRCA", 
#   data.category = "Gene expression",
#   data.type = "Gene expression quantification",
#   platform = "Illumina HiSeq", 
#   file.type  = "normalized_results",
#   experimental.strategy = "RNA-Seq",
#   legacy = TRUE
# )
# 
# GDCdownload(query)
# expr <- GDCprepare(query)
# 
# save(expr,file = "TCGA_BRCA_GDCprepare.RData")

# load("TCGA_BRCA_GDCprepare.RData")
# expr_dat <- assay(expr)


# ##TCGA BRCA clinical data
# query <- GDCquery(project = "TCGA-BRCA", 
#                   data.category = "Clinical",
#                   data.type = "Clinical Supplement", 
#                   data.format = "BCR Biotab")
# GDCdownload(query)
# clinical <- GDCprepare(query)
# save(clinical,file = "TCGA_BRCA_clinical_dat.RData")



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


#**************************************
#*no median centering
#*************************************
# #gene expr median centering
# tumor_intrinsic_expr_centering <- t(apply(tumor_intrinsic_expr,1,function(x){
#   x-median(x)
# }))


# #**********
# #*sigclust2
# #**********
# sigclust_res <- shc(tumor_intrinsic_expr_centering,metric="cor", linkage="average",null_alg = "2means")
# sigclust_res_2 <- shc(tumor_intrinsic_expr,metric="cor", linkage="average",null_alg = "2means")


#*********
#*sigclust
#*********

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



anno_color <- samples_of_classes %>% as.data.frame()
colnames(anno_color)[1] <- "Classes"
anno_color[,1] <- as.character(anno_color[,1])




#*************************************************************************
#*determine the correspondance between hclust label and intrinsic subtypes
#*************************************************************************
# #by heatmap?
# pheatmap(mat = expr_for_clust,
#          labels_row = rep("",nrow(expr_for_clust)),labels_col = rep("",ncol(expr_for_clust)),
#          angle_col = 45,color = colorRampPalette(c("green", "black", "firebrick3"))(101),
#          breaks=seq(-2,2,length.out=101),cluster_rows=gene_clust,cluster_cols = clust_res,
#          border_color=NA,annotation_col = anno_color)
# 
# pheatmap(mat = expr_for_clust[c("CDH3","CXCL1","KRT5","KRT17","CX3CL1","FZD7"),],
#          labels_row = rep("",nrow(expr_for_clust)),labels_col = rep("",ncol(expr_for_clust)),
#          angle_col = 45,color = colorRampPalette(c("green", "black", "firebrick3"))(101),
#          breaks=seq(-5,5,length.out=101),cluster_rows=FALSE,cluster_cols = clust_res,
#          border_color=NA,annotation_col = anno_color)



#***********************************************************************************
#****************Determine the cluster and subtype corresponding later**************
#****************Need to know which gene corresponds to which subtype***************
#***********************************************************************************
# #by counting expr?
# dat_for_corresponding <- expr_for_clust[intrinsic_gene_dat$intrinsic_genes,]
# 
# gene_rankings <- t(apply(X = expr_for_clust,1,function(x){
#   expr_rank <- rank(x,ties.method = "average")
#   tapply(expr_rank,INDEX = samples_of_classes,mean)
# }))
# 
# sample_gene_rankings <- t(apply(gene_rankings,1,rank))
# 
# cluster_gene_rankings <- apply(sample_gene_rankings,2,function(x){
#   ifelse(x==5,intrinsic_gene_dat$determined_subtype,"other")
# })
# 
# cluster_subtype_gene_proportion <- apply(apply(cluster_gene_rankings,2,table)[-6,],2,function(x){
#   x/table(intrinsic_gene_dat$determined_subtype)
# })
# 
# cluster_subtype_corresponding<- apply(cluster_subtype_gene_proportion,2,function(x){
#   names(which.max(x))
# })
# 
# #*****************
# #*for test scripts
# #*****************
# cluster_subtype_corresponding[1:5] <- c("LumA","LumB","Basal","Her2","Normal")
# if(length(unique(cluster_subtype_corresponding)) < 5)
#   stop("\nThe clusters in hierarchical clustering are not corresponded to five intrinsic subtypes!")
# 
# #cluster number to subtypes
# cluster_sample_subtypes <- as.factor(samples_of_classes)
# levels(cluster_sample_subtypes) <- cluster_subtype_corresponding
# cluster_sample_subtypes <- as.character(cluster_sample_subtypes)

cluster_sample_subtypes <- as.character(samples_of_classes)



#train model using AIMS idea
dat_for_model <- tumor_expr
datset_for_model <- list(x=dat_for_model,y=cluster_sample_subtypes)
datset_for_model_intrinsic_g <- list(x=tumor_intrinsic_expr,y=cluster_sample_subtypes)
res_prefix <- "TCGA_array_train_and_test"
res_prefix_intrinsic <- "TCGA_array_train_and_test_intrinsic_gene"

#using all genes
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Results/")
# train_and_test_with_label <- produce_train_test_set(expr_with_label = datset_for_model)
# train_and_test_res <- AIMS_train_and_pred(train_and_test_data_with_lable = train_and_test_with_label,
#                                           PREFIX = res_prefix,
#                                           k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res,file = paste0(res_prefix,"_AIMS_training_res.RData"))

# #load(paste0(res_prefix,"AIMS_training_res.RData"))

# #using intrinsic genes
train_and_test_with_label_intrinsic <- produce_train_test_set(expr_with_label = datset_for_model_intrinsic_g)
# train_and_test_res_intrinsic<-AIMS_train_and_pred(
#   train_and_test_data_with_lable=train_and_test_with_label_intrinsic,PREFIX = res_prefix_intrinsic,
#   k.fold = 10,num.of.rules = seq(1,50,1))
# save(train_and_test_res_intrinsic,file = paste0(res_prefix_intrinsic,"_AIMS_training_res.RData"))


#train model using PAM50 strategies
#using all genes
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
# #train model using PAM50 strategies
# #using intrinsic gens
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





































# #sigclust test for each cluster
# sig_res <- tapply(colnames(expr_for_clust),INDEX = samples_of_classes,function(x){
#   sigclust::sigclust(expr_for_clust[,x],nsim = 1000,icovest = 3)
#   #x
# },simplify = TRUE)
# sigclust_res_sig <- sigclust::sigclust(expr_for_clust,nsim = 1000,icovest = 3)
# 
# sigclust_res_sig <- sigclust::sigclust(expr_for_clust,nsim = 1000,icovest = 1)
# 
# 
# 
# 
# 
# plot(color_branches(clust_res,k=5),leaflab="none")
# plot(color_branches(clust_res,k=5),leaflab="textlike")
# plot(color_branches(clust_res,k=5),leaflab="perpendicular")



