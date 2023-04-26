#****************************************************************************
#*to compare two clustering strategy: pearson + average and euclidean+ ward.D
#****************************************************************************

library(TCGAbiolinks)
library(SummarizedExperiment)
#library(sigclust2)
library(tidyverse)
library(dendextend)
library(pheatmap)
library(viridis)
library(survival)
library(survminer)

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


# absolute expr to expr rank ----------------------------------------------

expr_for_clust <- apply(tumor_intrinsic_expr,2,function(x){
  rank(-1*x)
})







#do clust
#use pearson and average
cors <- cor(expr_for_clust,use = "pairwise.complete.obs", method = "pearson")
clust_res <- hclust(d = as.dist(1-cors),method = "average")

#cut to 5 classes
samples_of_classes <- cutree(clust_res,k = 5)
#plot(color_branches(clust_res,k=5),leaflab="none")

gene_clust <- hclust(d = dist(expr_for_clust,method = "euclidean"),method = "ward.D")

#do clust
clust_res_euclidean <- hclust(d = dist(t(expr_for_clust),method = "euclidean"),method = "ward.D")
samples_of_classes_euclidean <- cutree(clust_res_euclidean,k=5)

#do clust by median centering
clust_dat_centering <- t(apply(expr_for_clust,1,function(x){
  x-median(x)
}))

cors_centering <- cor(clust_dat_centering,use = "pairwise.complete.obs", method = "pearson")
clust_res_centering <- hclust(d = as.dist(1-cors_centering),method = "average")
samples_of_classes_centering <- cutree(clust_res_centering,k=5)





cluster_biomarker_res <- samples_of_classes_euclidean %>% as.data.frame() %>% 
  mutate(barcode=names(samples_of_classes_euclidean)) %>% rename(cluster_euclidean=1) %>% 
  left_join(A,by = "barcode") %>% 
  left_join(samples_of_classes %>% as.data.frame() %>%
              mutate(barcode=names(samples_of_classes)) %>% rename(cluster_pearson=1),by="barcode") %>%
  left_join(samples_of_classes_centering %>% as.data.frame() %>%
              mutate(barcode=names(samples_of_classes_centering)) %>% 
              rename(cluster_pearson_centering=1),by="barcode")

cluster_biomarker_res <- cluster_biomarker_res %>% 
  select(barcode,patient,cluster_pearson,cluster_euclidean,cluster_pearson_centering,ajcc_pathologic_stage,
         days_to_last_follow_up,primary_diagnosis,paper_BRCA_Subtype_PAM50,ends_with("Clusters"))





load("TCGA_BRCA_clinical_dat.RData")
#clinical$clinical_patient_brca %>% full_join(cluster_biomarker_res,by = c("bcr_patient_barcode"="patient"))


cluster_biomarker_res %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>%
  select(barcode,patient,cluster_pearson,cluster_euclidean,cluster_pearson_centering,paper_BRCA_Subtype_PAM50,
         ajcc_pathologic_stage,days_to_last_follow_up,primary_diagnosis,er_status_by_ihc,
         pr_status_by_ihc,her2_status_by_ihc,vital_status,ends_with("Clusters")) %>% 
  arrange(cluster_euclidean,cluster_pearson,cluster_pearson_centering)

dat_of_clinical <- cluster_biomarker_res %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>%
  select(barcode,patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson,paper_BRCA_Subtype_PAM50,
         ajcc_pathologic_stage,days_to_last_follow_up,primary_diagnosis,er_status_by_ihc,
         pr_status_by_ihc,her2_status_by_ihc,vital_status,ends_with("Clusters")) %>% 
  arrange(cluster_euclidean,cluster_pearson,cluster_pearson_centering) %>% 
  apply(2,as.character) %>% as.data.frame()

# dat_of_clinical %>%
#   select(patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson,paper_BRCA_Subtype_PAM50,
#          er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc)  %>%
#   pivot_longer(cols = 2:last_col(),
#                names_to = "metrics",values_to = "value") %>%
#   mutate(value=ifelse(grepl("cluster",metrics),paste0(metrics,"_",value),value)) %>%
#   mutate(value=ifelse(grepl("Not",value),NA,value)) %>%
#   mutate(value = factor(value,
#                         levels = c(paste0(rep("cluster_",15),
#                                           rep(c("euclidean","pearson","pearson_centering"),each=5),"_",1:5),
#                                    c("Basal","Her2","LumA","LumB","Normal","Positive","Negative",
#                                      "Indeterminate","Equivocal")))) %>%
#   mutate(value = as.integer(value)) %>%
#   arrange(metrics) %>% as.data.frame() %>%
#   ggplot() +
#   geom_tile(aes(x=patient,y=metrics,fill=value)) +
#   scale_fill_viridis(discrete=FALSE) +
#   theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
#   labs(x='',y='',fill='Coef in ARD')



#compare 3 clustering res
dat_of_clinical %>% 
  select(patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson)  %>%
  pivot_longer(cols = 2:last_col(),
               names_to = "metrics",values_to = "value") %>% 
  mutate(value=ifelse(grepl("cluster",metrics),paste0(metrics,"_",value),value)) %>%
  mutate(value=ifelse(grepl("Not",value),NA,value)) %>% 
  mutate(value = factor(value,levels = c(paste0(rep("cluster_",15),
                                                rep(c("euclidean","pearson","pearson_centering"),each=5),"_",1:5),
                                         c("Basal","Her2","LumA","LumB","Normal","Positive","Negative",
                                           "Indeterminate","Equivocal")))) %>%
  mutate(value = as.integer(value)) %>%
  arrange(value) %>% mutate(patient = factor(patient,levels = unique(patient))) %>% 
  ggplot() +
  geom_tile(aes(x=patient,y=metrics,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Set3"),"red","green","blue")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
  labs(x='',y='',fill='Coef in ARD')

#compare euclidean and pam50 ER PR Her2
dat_of_clinical %>% 
  select(patient,cluster_euclidean,paper_BRCA_Subtype_PAM50,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc) %>%
  pivot_longer(cols = 2:last_col(),
               names_to = "metrics",values_to = "value") %>%
  mutate(value=ifelse(grepl("cluster",metrics),paste0(metrics,"_",value),value)) %>%
  mutate(value=ifelse(grepl("Not",value),NA,value)) %>%
  mutate(value = factor(value,levels = c(paste0(rep("cluster_",15),
                                                rep(c("euclidean","pearson","pearson_centering"),each=5),"_",1:5),
                                         c("Basal","Her2","LumA","LumB","Normal","Positive","Negative"))))%>%
  mutate(value = as.integer(value)) %>%
  arrange(value) %>% mutate(patient = factor(patient,levels = unique(patient)))%>% 
  mutate(metrics=factor(metrics,levels = c("cluster_euclidean","cluster_pearson_centering","cluster_pearson",
                                           "paper_BRCA_Subtype_PAM50","er_status_by_ihc","pr_status_by_ihc",
                                           "her2_status_by_ihc"))) %>% 
  ggplot() +
  geom_tile(aes(x=patient,y=metrics,fill=factor(value,levels = unique(value)))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(10,"Set3"),"black","white")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  labs(x='',y='',fill='Coef in ARD')


#compare pearson centering and pam50 ER PR Her2
dat_of_clinical %>% 
  select(patient,cluster_pearson_centering,paper_BRCA_Subtype_PAM50,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc) %>%
  pivot_longer(cols = 2:last_col(),
               names_to = "metrics",values_to = "value") %>%
  mutate(value=ifelse(grepl("cluster",metrics),paste0(metrics,"_",value),value)) %>%
  mutate(value=ifelse(grepl("Not",value),NA,value)) %>%
  mutate(value = factor(value,levels = c(paste0(rep("cluster_",15),
                                                rep(c("euclidean","pearson","pearson_centering"),each=5),"_",1:5),
                                         c("Basal","Her2","LumA","LumB","Normal","Positive","Negative")))) %>%
  mutate(value = as.integer(value)) %>%
  arrange(value) %>% mutate(patient = factor(patient,levels = unique(patient)))%>% 
  mutate(metrics=factor(metrics,levels = c("cluster_euclidean","cluster_pearson_centering","cluster_pearson",
                                           "paper_BRCA_Subtype_PAM50","er_status_by_ihc","pr_status_by_ihc",
                                           "her2_status_by_ihc"))) %>% 
  ggplot() +
  geom_tile(aes(x=patient,y=metrics,fill=factor(value,levels = unique(value)))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(10,"Set3"),"black","white")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  labs(x='',y='',fill='Coef in ARD')


#compare pearson and pam50 ER PR Her2
dat_of_clinical %>% 
  select(patient,cluster_pearson,paper_BRCA_Subtype_PAM50,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc) %>%
  pivot_longer(cols = 2:last_col(),
               names_to = "metrics",values_to = "value") %>%
  mutate(value=ifelse(grepl("cluster",metrics),paste0(metrics,"_",value),value)) %>%
  mutate(value=ifelse(grepl("Not",value),NA,value)) %>%
  mutate(value = factor(value,levels = c(paste0(rep("cluster_",15),
                                                rep(c("euclidean","pearson","pearson_centering"),each=5),"_",1:5),
                                         c("Basal","Her2","LumA","LumB","Normal","Positive","Negative")))) %>%
  mutate(value = as.integer(value)) %>%
  arrange(value) %>% mutate(patient = factor(patient,levels = unique(patient))) %>% 
  mutate(metrics=factor(metrics,levels = c("cluster_euclidean","cluster_pearson_centering","cluster_pearson",
                                           "paper_BRCA_Subtype_PAM50","er_status_by_ihc","pr_status_by_ihc",
                                           "her2_status_by_ihc"))) %>%
  ggplot() +
  geom_tile(aes(x=patient,y=metrics,fill=factor(value,levels = unique(value)))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(10,"Set3"),"black","white")) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  labs(x='',y='',fill='Coef in ARD')


##survival analysis
#euclidean clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~cluster_euclidean,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#pearson centering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~cluster_pearson_centering,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#pearson centering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~cluster_pearson,data = surv_dat)
ggsurvplot(fit,pval = TRUE)



#subtyping/clustering res in Nature
#paper PAM50 subtyping
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_BRCA_Subtype_PAM50,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#mRNA clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_mRNA.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#CNV clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_CNV.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#Mutation clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_Mutation.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#DNA Methylation clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_DNA.Methylation.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#miRNA clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_miRNA.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#lncRNA clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_lncRNA.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#protein clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_Protein.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#PARADIGM clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_PARADIGM.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#Pan.Gyn clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",1,0))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_Pan.Gyn.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)



#sigclust_res_sig <- sigclust::sigclust(expr_for_clust,nsim = 1000,icovest = 1)





#Nature 2012 Breast Cancer data
#Clinical data
Nature_2012_clinical <- read.table("2012_Nature_BRCA_microarray_clinical_data.txt",header = TRUE,sep = "\t",
                                   stringsAsFactors = FALSE)
Nature_2012_clinical %>% rename(patient="Complete.TCGA.ID") %>% select(patient,PAM50.mRNA) %>% 
  filter(!is.na(PAM50.mRNA))

Nature_2012_clinical %>% rename(patient="Complete.TCGA.ID") %>% select(patient,PAM50.mRNA) %>%
  filter(!is.na(PAM50.mRNA)) %>% 
  left_join(dat_of_clinical,by = "patient")

a <- Nature_2012_clinical %>% rename(patient="Complete.TCGA.ID") %>% select(patient,PAM50.mRNA) %>%
  filter(!is.na(PAM50.mRNA)) %>% 
  left_join(dat_of_clinical,by = "patient")

#survival
surv_dat <- Nature_2012_clinical %>% select(starts_with("OS"),PAM50.mRNA)
fit <- survfit(Surv(OS.Time,OS.event)~PAM50.mRNA,data = surv_dat)
ggsurvplot(fit,pval = TRUE)



##Using new subtype to subtype TCGA samples

#load model
load("../Results/TCGA_train_and_testAIMS_training_res.RData")


setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
source("clanc.R")
source("../../New_PAM50_classifier/scripts/all_functions_20220302.R")
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/New_PAM50_classifier/scripts/trainAIMS-master/trainAIMS-master/")
source("trainAIMS_2.R")

setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")
a <- predict.one.vs.all.tsp(D = tumor_expr,GeneName = rownames(tumor_expr),
                            one.vs.all.tsp = AIMS_training_res$training_res$final_model)











