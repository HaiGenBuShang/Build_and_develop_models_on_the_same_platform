library(TCGAbiolinks)
library(SummarizedExperiment)
library(sigclust2)
library(tidyverse)
library(dendextend)
library(pheatmap)
library(pamr)
library(seqc)


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

##function to select genes
select_genes <- function(expr_dt,expr_pct=0.7,expr_thres=0){
  n_clmn <- ncol(expr_dt)
  #if the expression of one gene in a samples is bigger than expr_thres
  total_count <- apply(expr_dt,1,function(x){
    sum(x > expr_thres)
  })
  total_count[rownames(expr_dt)] > expr_pct*n_clmn
}







setwd("../Data/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/")
setwd("./Gene_Expression_Quantification/0019c951-16c5-48d0-85c8-58d96b12d330/")
Ensembl_and_Entrez <- read.table("ba295155-272e-43eb-9d6a-e4c9c392e68b.rna_seq.augmented_star_gene_counts.tsv",
                                 header = TRUE,
                                 sep = "\t",stringsAsFactors = FALSE)[-1:-4,1:3]
Ensembl_and_Entrez$gene_id <- str_remove_all(Ensembl_and_Entrez$gene_id,"\\..*")



ILM_refseq_gene_AGR <- ILM_refseq_gene_AGR[!is.na(ILM_refseq_gene_AGR$Symbol),]

SEQC_ILM_AGR_FPKM_entrez <- apply(ILM_refseq_gene_AGR[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_AGR[,3])

SEQC_ILM_AGR_FPKM_symbol <- apply(SEQC_ILM_AGR_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_AGR[,2])

SEQC_expr_AGR <- log2(SEQC_ILM_AGR_FPKM_symbol[intersect(Ensembl_and_Entrez$gene_name,
                                                         rownames(SEQC_ILM_AGR_FPKM_symbol)),]+0.01)


ILM_refseq_gene_BGI <- ILM_refseq_gene_BGI[!is.na(ILM_refseq_gene_BGI$Symbol),]
SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_BGI[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_BGI[,3])


SEQC_ILM_BGI_FPKM_symbol <- apply(SEQC_ILM_BGI_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_BGI[,2])


SEQC_expr_BGI <- log2(SEQC_ILM_BGI_FPKM_symbol[intersect(Ensembl_and_Entrez$gene_name,
                                                         rownames(SEQC_ILM_BGI_FPKM_symbol)),]+0.01)






# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("../../../../../../../Data/")

#tumor sample expr
tumor_expr <- SEQC_expr_AGR
tumor_expr <- tumor_expr[apply(tumor_expr,1,function(x)sd(x)!=0),]

#intrinsic genes
intrinsic_genes <- read.table("2009_JCO_intrinsic_genes_S_table_5.txt",header = FALSE,sep = "\t",
                              stringsAsFactors = FALSE)
#tumor intrinsic gene exprs
#****
#*use the direct intersect between intrinsic genes and genes in dataset,
#*which need to be done more rigorously
#****
tumor_intrinsic_expr <- tumor_expr[intersect(intrinsic_genes[,1],rownames(tumor_expr)),]


RNA_seq_expr_for_clust <- tumor_intrinsic_expr
RNA_seq_expr <- tumor_expr


#tumor sample expr
tumor_expr <- SEQC_expr_BGI
tumor_expr <- tumor_expr[apply(tumor_expr,1,function(x)sd(x)!=0),]

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
RNA_seq_expr_PAM50 <- RNA_seq_expr[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]


intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,scale))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,scale))

combined_dat <- cbind(RNA_seq_expr_PAM50[intersect_genes,],array_expr_PAM50[intersect_genes,])
combined_dat <- combined_dat[!apply(combined_dat,1,function(x){is.na(x) %>% any()}),]


colnames(combined_dat) <- c(paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_AGR"),
                            paste0(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                   "_BGI"))


#clustering methods were changed from ward.D to Ward.D
clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D")
#gene_clustering_res <- hclust(d = dist(combined_dat),method = "ward.D")



# samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
samples_of_classes <- cutree(clustering_res,
                             k = combined_dat %>% colnames() %>% str_remove_all("_.*") %>% unique() %>% length()) %>% 
  as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1) %>% mutate(cluster_num=as.character(cluster_num))


consensus_samples_array_and_RNA_seq <- samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>%
  group_by(samples,cluster_num) %>% count() %>% filter(n==128)

# samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
#   group_by(samples,cluster_num) %>% count() %>% filter(n!=2)





samples_of_5_classes <- cutree(clustering_res,
                               k=combined_dat %>% colnames() %>% str_remove_all("_.*") %>% unique() %>% length()) %>% 
  as.matrix() %>% as.data.frame() %>% 
  rownames_to_column("barcode") %>% rename(cluster_num=2) %>% mutate(cluster_num=as.character(cluster_num)) %>% 
  mutate(samples=str_replace_all(barcode,"_.*",""))

consensus_samples_of_5_classes <- samples_of_5_classes %>% 
  inner_join(consensus_samples_array_and_RNA_seq,by="samples") %>% 
  mutate(data_type=str_replace_all(barcode,".*_","")) %>% mutate(sample_name=str_remove_all(barcode,"_AGR|_BGI")) %>% 
  pivot_wider(id_cols = -barcode,names_from = data_type,values_from = data_type)



# train model using AIMS idea ---------------------------------------------

#intersect genes
intersect_intrinsic_genes <- intersect(intersect_genes,intrinsic_genes[,1])
intersect_all_genes <- intersect(rownames(RNA_seq_expr),rownames(array_expr))





#*************
#*RNA-seq data
#*************
tumor_expr <- RNA_seq_expr[intersect_all_genes,consensus_samples_of_5_classes$sample_name]

#****drop genes with the same expression in at least 30% samples
#****other wise the clanc training might cause error
tumor_expr <- tumor_expr[apply(tumor_expr,1,function(x){
  (table(x) %>% max)<=(ncol(tumor_expr)*0.30)
}),]



cluster_sample_subtypes <- consensus_samples_of_5_classes$cluster_num.x
tumor_intrinsic_expr <- RNA_seq_expr[intersect_intrinsic_genes,consensus_samples_of_5_classes$sample_name]
#****drop genes with the same expression in at least 30% samples
#****other wise the clanc training might cause error
tumor_intrinsic_expr <- tumor_intrinsic_expr[apply(tumor_intrinsic_expr,1,function(x){
  (table(x) %>% max)<=(ncol(tumor_intrinsic_expr)*0.30)
}),]




dat_for_model <- tumor_expr
datset_for_model <- list(x=dat_for_model,y=cluster_sample_subtypes)
datset_for_model_intrinsic_g <- list(x=tumor_intrinsic_expr,y=cluster_sample_subtypes)
save(datset_for_model,datset_for_model_intrinsic_g,
     file = "../Results/SEQC_datset_for_model_and_dataset_for_model_intrinsic_g.RData")

res_prefix <- "SEQC_AGR_BGI_Z_score_AGR_model"
res_prefix_intrinsic <- "SEQC_AGR_BGI_Z_score_AGR_model_intrinsic_gene"

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
# #using all genes
datset_for_model$x <- tumor_expr[sample(1:nrow(tumor_expr),2000),]
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
# train model using PAM50 strategies
# using all genes
set.seed(12345678)
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















