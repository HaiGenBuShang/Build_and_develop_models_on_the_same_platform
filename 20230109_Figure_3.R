#****************************************************************************
#*all figures in manuscript by Fig. 3
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



# Figure 3 ----------------------------------------------------------------

rm(list = ls())

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

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

# pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),labels_col = rep('',ncol(combined_dat)),
#          cutree_cols = 5,fontsize_col = 2,clustering_method = "ward.D",
#          color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101),
#          cluster_cols = TRUE,scale = "row",cluster_rows = TRUE)

#discrimination genes in PNAS 2003 Sorlie Figure 1
Sorlie_genes <- c("TLK1","TRAP100","PPARBP","ERBB2","GRB7","ERBB2","ATP5G1","PRNPIP","NSEP1","GGH","LAPTM4B","PRDX4",
                  "CCNE1","SQLE","CXCL1","CDH3","ANXA8","KRT5","TRIM29","KRT17","MFGE8","CX3CL1","FZD7","CHI3L2",
                  "B3GNT5",
                  "PIK3R1","AKR1C1","FACL2","LRBA","NAT1","LIV-1","HNF3A","XBP1","GATA3","ESR1","PTP4A2","RERG",
                  "SCUBE2") %>% unique()

#clustering methods were changed from ward.D to Ward.D
clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D")
gene_clustering_res <- hclust(d = dist(combined_dat),method = "ward.D")


pdf("../Results/Figure/Figure_3_supp_1.pdf",height = 50,width = 50)
pheatmap(mat = combined_dat,
         labels_row = ifelse(rownames(combined_dat)%in%Sorlie_genes,rownames(combined_dat),""),
         labels_col = rep('',ncol(combined_dat)),cutree_cols = 5,
         fontsize_col = 2,clustering_method = "ward.D",
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101),
         cluster_cols = clustering_res,scale = "row",cluster_rows = TRUE)
dev.off()


data("pam50")
RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(RNA_seq_expr_PAM50),
                                           annot = data.frame(Gene.Symbol=rownames(RNA_seq_expr_PAM50),
                                                              stringsAsFactors = FALSE),
                                           do.mapping = FALSE)
array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(array_expr_PAM50),
                                         annot = data.frame(Gene.Symbol=rownames(array_expr_PAM50),
                                                            stringsAsFactors = FALSE),
                                         do.mapping = FALSE)

anno_cols <- data.frame(PAM50=c(RNA_seq_PAM50$subtype,array_PAM50$subtype),
                        row.names = c(paste0(intersect(colnames(RNA_seq_expr_for_clust),
                                                       colnames(array_expr_for_clust)),
                                             "_RNAseq"),
                                      paste0(intersect(colnames(RNA_seq_expr_for_clust),
                                                       colnames(array_expr_for_clust)),
                                             "_array")))

pdf("../Results/Figure/Figure_3_supp_2.pdf",height = 10,width = 10)
pheatmap(mat = combined_dat[intersect(rownames(combined_dat),Sorlie_genes),],
         labels_row = intersect(rownames(combined_dat),Sorlie_genes),
         labels_col = rep('',ncol(combined_dat)),cutree_cols = 5,cutree_rows = 5,
         fontsize_col = 2,clustering_method = "ward.D",
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101),
         cluster_cols = clustering_res,scale = "row",cluster_rows = TRUE,annotation_col = anno_cols)
dev.off()

samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1) %>% mutate(cluster_num=as.character(cluster_num))


consensus_samples_array_and_RNA_seq <- samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n!=2)


(samples_of_5_classes <- cutree(clustering_res,k=5) %>% as.matrix() %>% as.data.frame() %>% 
  rownames_to_column("barcode") %>% rename(cluster_num=2) %>% mutate(cluster_num=as.character(cluster_num)) %>% 
  mutate(samples=str_replace_all(barcode,"_.*","")))

consensus_samples_of_5_classes <- samples_of_5_classes %>% 
  inner_join(consensus_samples_array_and_RNA_seq,by="samples") %>% 
  mutate(data_type=str_replace_all(barcode,".*_","")) %>% 
  pivot_wider(id_cols = -barcode,names_from = data_type,values_from = data_type)



# #sample clustering order
# #******************************************
# clustering_res$labels[clustering_res$order]
# #******************************************

#clinical data
load("TCGA_BRCA_clinical_dat.RData")

#heatmap of RNA_seq and expression

RNA_seq_dat_for_clustering <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),
                                                                colnames(array_expr_for_clust))]
RNA_samples <- str_replace_all(grep("_RNAseq",clustering_res$labels[clustering_res$order],value = TRUE),"_RNA.*","")
combined_dat <- RNA_seq_dat_for_clustering[gene_clustering_res$labels[gene_clustering_res$order],
                                           RNA_samples[RNA_samples%in%consensus_samples_of_5_classes$samples]]

annotation_dat <- consensus_samples_of_5_classes %>% select(1,2) %>% rename(Subtype=1) %>% 
  rename(Subtype=1,barcode=2) %>% left_join(A,by = "barcode") %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>% 
  select(Subtype,barcode,er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc) %>% 
  rename(ER_status=er_status_by_ihc,PR_status=pr_status_by_ihc,HER2_status=her2_status_by_ihc) %>% 
  mutate(ER_status=ifelse(grepl("Positive|Negative",ER_status),ER_status,NA)) %>% 
  mutate(PR_status=ifelse(grepl("Positive|Negative",PR_status),PR_status,NA)) %>% 
  mutate(HER2_status=ifelse(grepl("Positive|Negative",HER2_status),HER2_status,NA))

annotation_col <- data.frame(annotation_dat,row.names = 2)
anno_color <- list(ER_status=c(Positive="black",Negative = "grey"),
                   PR_status=c(Positive="black",Negative = "grey"),
                   HER2_status=c(Positive="black",Negative = "grey"),
                   Subtype=set_names(RColorBrewer::brewer.pal(5,"Set1"),as.character(1:5)))

pdf("../Results/Figure/Figure_3_RNA_seq_expr_and_anno.pdf",height = 5,width = 10)
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),labels_col = rep('',ncol(combined_dat)),
         annotation_col = annotation_col,annotation_colors = anno_color,fontsize_col = 2,clustering_method = "ward.D",
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101),
         cluster_cols = FALSE,scale = "row",cluster_rows = FALSE)
dev.off()

array_for_clustering <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),
                                                        colnames(array_expr_for_clust))]
array_samples <- str_replace_all(grep("_array",clustering_res$labels[clustering_res$order],value = TRUE),"_array.*","")
combined_dat <- array_for_clustering[gene_clustering_res$labels[gene_clustering_res$order],
                                     array_samples[array_samples%in%consensus_samples_of_5_classes$samples]]

annotation_dat <-consensus_samples_of_5_classes %>% select(1,2) %>% rename(Subtype=1) %>% 
  rename(Subtype=1,barcode=2) %>% left_join(A,by = "barcode") %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>% 
  select(Subtype,barcode,er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc) %>% 
  rename(ER_status=er_status_by_ihc,PR_status=pr_status_by_ihc,HER2_status=her2_status_by_ihc) %>% 
  mutate(ER_status=ifelse(grepl("Positive|Negative",ER_status),ER_status,NA)) %>% 
  mutate(PR_status=ifelse(grepl("Positive|Negative",PR_status),PR_status,NA)) %>% 
  mutate(HER2_status=ifelse(grepl("Positive|Negative",HER2_status),HER2_status,NA))

annotation_col <- data.frame(annotation_dat,row.names = 2)
anno_color <- list(ER_status=c(Positive="black",Negative = "grey"),
                   PR_status=c(Positive="black",Negative = "grey"),
                   HER2_status=c(Positive="black",Negative = "grey"),
                   Subtype=set_names(RColorBrewer::brewer.pal(5,"Set1"),as.character(1:5)))

pdf("../Results/Figure/Figure_3_array_expr_and_anno.pdf",height = 5,width = 10)
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),labels_col = rep('',ncol(combined_dat)),
         annotation_col = annotation_col,annotation_colors = anno_color,fontsize_col = 2,clustering_method = "ward.D",
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101),
         cluster_cols = FALSE,scale = "row",cluster_rows = FALSE)
dev.off()



# BRCA prototype oncoplot ----------------------------------------------------------------

#
#clinical data
load("TCGA_BRCA_clinical_dat.RData")

cluster_biomarker_res <- consensus_samples_of_5_classes %>% select(1:2) %>% rename(Subtype=1,barcode=2) %>%
  bind_rows(data.frame(Subtype=NA,
                      barcode=setdiff(intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust)),
                                      consensus_samples_of_5_classes$samples))) %>%
  left_join(A,by = "barcode")

cluster_biomarker_res <- cluster_biomarker_res %>%
  select(Subtype,barcode,patient,ajcc_pathologic_stage,
         days_to_last_follow_up,primary_diagnosis,paper_BRCA_Subtype_PAM50,ends_with("Clusters"))

dat_of_clinical <- cluster_biomarker_res %>%
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>%
  select(Subtype,barcode,patient,paper_BRCA_Subtype_PAM50,
         ajcc_pathologic_stage,days_to_last_follow_up,death_days_to,primary_diagnosis,er_status_by_ihc,
         pr_status_by_ihc,her2_status_by_ihc,vital_status,ends_with("Clusters")) %>%
  apply(2,as.character) %>% as.data.frame() %>%
  mutate(days_to_death=str_replace_all(death_days_to,'[^0-9]','')) %>%
  mutate(days_to_last_follow_up=ifelse(days_to_death=='',days_to_last_follow_up,days_to_death)) %>%
  mutate(days_to_death=ifelse(days_to_death=='',NA,as.integer(days_to_death))) %>%
  mutate(days_to_last_follow_up=as.integer(days_to_last_follow_up)) %>% filter(!is.na(vital_status))

dat_of_clinical <- dat_of_clinical %>%
  rename(ER_status=er_status_by_ihc,PR_status=pr_status_by_ihc,HER2_status=her2_status_by_ihc) %>%
  mutate(ER_status=ifelse(grepl("Positive|Negative",ER_status),ER_status,NA)) %>%
  mutate(PR_status=ifelse(grepl("Positive|Negative",PR_status),PR_status,NA)) %>%
  mutate(HER2_status=ifelse(grepl("Positive|Negative",HER2_status),HER2_status,NA))


#mutation info
load("20221215_BRCA_mutation.RData")

#oncoplot by maftools
maftools_maf <- read.maf(maf = maf %>%
                           mutate(Tumor_Sample_Barcode=str_replace_all(Tumor_Sample_Barcode,
                                                                       "(TCGA)-(..)-(....).*","\\1-\\2-\\3")) %>%
                           filter(Tumor_Sample_Barcode%in%dat_of_clinical$patient),
                         clinicalData = dat_of_clinical %>% rename(Tumor_Sample_Barcode=patient))
pdf("../Results/Figure/Figure_3_oncoplot.pdf",width = 10,height = 7.5)
oncoplot(
  maf = maftools_maf, #genes = c("BRCA1","BRCA2"),
  clinicalFeatures = c('Subtype','ER_status','PR_status',"HER2_status"),
  top = 20,
  sortByAnnotation = TRUE,removeNonMutated = FALSE,
  annotationColor = list(ER_status=c(Positive="black",Negative = "grey"),
                         PR_status=c(Positive="black",Negative = "grey"),
                         HER2_status=c(Positive="black",Negative = "grey"),
                         Subtype=set_names(RColorBrewer::brewer.pal(5,"Set1"),as.character(1:5))),
  fontSize = 0.8,anno_height = 2)
dev.off()


dat_of_clinical %>% group_by(Subtype,ER_status,PR_status,HER2_status) %>% count() %>% group_by(Subtype) %>%
  mutate(total=sum(n)) %>% mutate(triN_rate=n/total) %>% arrange(desc(triN_rate))




#Define immunohistology subtype
imhi_subtype_def <- list(LumA=c("Positive, Positive, Negative",
                                "Positive, [Not Evaluated], Negative",
                                "[Not Evaluated], Positive, Negative",
                                "Positive, Negative, Negative",
                                "Negative, Positive, Negative"),
                         LumB=c("Positive, Positive, Positive",
                                "Positive, [Not Evaluated], Positive",
                                "[Not Evaluated], Positive, Positive",
                                "Positive, Negative, Positive",
                                "Negative, Positive, Positive"),
                         Her2=c("Negative, Negative, Positive"),
                         Basal=c(c("Negative, Negative, Negative")))

eval_imhi_subtype <- function(imhi_info,imhi_subtype_def){
  res <- names(imhi_subtype_def)[map_lgl(imhi_subtype_def,~(any(imhi_info%in%.)))]
  ifelse(length(res)<1,"Other",res)
}

annotation_col %>% rowwise() %>% mutate(imhi=paste0(c_across(2:4),collapse = ", ")) %>% 
  mutate(imhi_subtype = map_chr(imhi,eval_imhi_subtype,imhi_subtype_def)) %>% group_by(Subtype) %>% count(imhi_subtype)



#Survival analysis
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient)) %>% 
  rename(PAM50=paper_BRCA_Subtype_PAM50)
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~PAM50,data = surv_dat)
g_surv <- ggsurvplot(fit,pval = TRUE,risk.table = TRUE,fontsize=7,
                     font.tickslab = 15)
# file.remove(filename = "../Results/Figure/Figure_3_survival_PAM50.pdf")
ggsave(filename = "../Results/Figure/Figure_3_survival_PAM50.pdf",plot = g_surv,
       width = 10,height = 10,device = "pdf")



surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient)) %>% 
  rename(PAM50=paper_BRCA_Subtype_PAM50)
fit2 <- survfit(Surv(days_to_last_follow_up,vital_status)~Subtype,data = surv_dat)
g_surv <- ggsurvplot(fit2,pval = TRUE,risk.table = TRUE,fontsize=7,
                     font.tickslab = 15)
ggsave(filename = "../Results/Figure/Figure_3_survival_Subtype.pdf",plot = g_surv,width = 10,height = 10,
       device = "pdf")





















