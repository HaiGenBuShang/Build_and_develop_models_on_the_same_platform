#****************************************************************************
#*to compare two clustering strategy: pearson + average and euclidean+ ward.D
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
expr_for_clust <- tumor_intrinsic_expr

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
         ajcc_pathologic_stage,days_to_last_follow_up,death_days_to,primary_diagnosis,er_status_by_ihc,
         pr_status_by_ihc,her2_status_by_ihc,vital_status,ends_with("Clusters")) %>% 
  arrange(cluster_euclidean,cluster_pearson,cluster_pearson_centering)

dat_of_clinical <- cluster_biomarker_res %>% 
  left_join(clinical$clinical_patient_brca,by = c("patient"="bcr_patient_barcode")) %>%
  select(barcode,patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson,paper_BRCA_Subtype_PAM50,
         ajcc_pathologic_stage,days_to_last_follow_up,death_days_to,primary_diagnosis,er_status_by_ihc,
         pr_status_by_ihc,her2_status_by_ihc,vital_status,ends_with("Clusters")) %>% 
  arrange(cluster_euclidean,cluster_pearson,cluster_pearson_centering) %>% 
  apply(2,as.character) %>% as.data.frame() %>% 
  mutate(days_to_death=str_replace_all(death_days_to,'[^0-9]','')) %>%
  mutate(days_to_last_follow_up=ifelse(days_to_death=='',days_to_last_follow_up,days_to_death)) %>%
  mutate(days_to_death=ifelse(days_to_death=='',NA,as.integer(days_to_death))) %>%
  mutate(days_to_last_follow_up=as.integer(days_to_last_follow_up)) %>% filter(!is.na(vital_status))
  
##Pam50 subtype of the same patient in Nature paper
dat_of_clinical %>% group_by(patient) %>% 
  nest() %>%  
  mutate(number = map_dbl(data,nrow)) %>% filter(number>1) %>% 
  mutate(subtype = map_chr(data,~(paste0(.$paper_BRCA_Subtype_PAM50,collapse = ", "))))

##subtype of the same patient in cluster euclidean
dat_of_clinical %>% group_by(patient) %>% 
  nest() %>%  
  mutate(number = map_dbl(data,nrow)) %>% filter(number>1) %>% 
  mutate(subtype = map_chr(data,~(paste0(.$cluster_euclidean,collapse = ", "))))

##subtype of the same patient in cluster pearson centering
dat_of_clinical %>% group_by(patient) %>% 
  nest() %>%  
  mutate(number = map_dbl(data,nrow)) %>% filter(number>1) %>% 
  mutate(subtype = map_chr(data,~(paste0(.$cluster_pearson_centering,collapse = ", "))))

##subtype of the same patient in cluster pearson centering
dat_of_clinical %>% group_by(patient) %>% 
  nest() %>%  
  mutate(number = map_dbl(data,nrow)) %>% filter(number>1) %>% 
  mutate(subtype = map_chr(data,~(paste0(.$cluster_pearson,collapse = ", "))))




#compare 3 clustering res
dat_of_clinical %>% 
  select(patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson) %>%
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



dat_of_clinical %>% 
  select(patient,cluster_euclidean,paper_BRCA_Subtype_PAM50,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc) %>% group_by(cluster_euclidean) %>% 
  count(paper_BRCA_Subtype_PAM50)










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
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~cluster_euclidean,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE)
#pearson centering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~cluster_pearson_centering,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#pearson centering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~cluster_pearson,data = surv_dat)
ggsurvplot(fit,pval = TRUE)



#subtyping/clustering res in Nature
#paper PAM50 subtyping
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_BRCA_Subtype_PAM50,data = surv_dat)
ggsurvplot(fit,pval = TRUE,risk.table = TRUE)
#mRNA clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_mRNA.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#CNV clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_CNV.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#Mutation clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_Mutation.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#DNA Methylation clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_DNA.Methylation.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#miRNA clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_miRNA.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#lncRNA clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_lncRNA.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#protein clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_Protein.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#PARADIGM clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
fit <- survfit(Surv(days_to_last_follow_up,vital_status)~paper_PARADIGM.Clusters,data = surv_dat)
ggsurvplot(fit,pval = TRUE)
#Pan.Gyn clustering
surv_dat <- dat_of_clinical %>% mutate(days_to_last_follow_up=as.numeric(days_to_last_follow_up)/365) %>%
  mutate(vital_status=ifelse(vital_status=="Alive",0,1)) %>% filter(!duplicated(patient))
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






# TCGA PAM50 subtying by myself -------------------------------------------

load("/mnt/Miscrosoft/TCGA_project/Data/R_data/TCGA-BRCA/Expr_and_pheno.RData")
expr_dat <- X1_genename
expr_dat <- log2(expr_dat+1)
expr_dat <- expr_dat[!duplicated(rownames(expr_dat)),]
tumor_sample_ID_dat <- A %>% 
  dplyr::filter(sample_type=="Primary Tumor") %>% dplyr::select(barcode,sample_type)

#tumor sample expr
tumor_expr <- expr_dat[,tumor_sample_ID_dat$barcode]
data("pam50.robust")
#note that different normalization method would affect the results
PAM50_subtyping_res <- molecular.subtyping(sbt.model = "pam50",data = t(tumor_expr),
                                           annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                              stringsAsFactors = FALSE),
                                           do.mapping = FALSE)

new_dat_of_clinidal <- dat_of_clinical %>% 
  left_join(as.data.frame(PAM50_subtyping_res$subtype) %>% mutate(barcode=rownames(.)) %>% 
              rename(subtype_PAM50_classifier=1), by = "barcode") %>% 
  select(paper_BRCA_Subtype_PAM50,subtype_PAM50_classifier,
         patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc,barcode)




# New prototype characterisctic -------------------------------------------

#the different samples from the same patient are in the sample Immunohist status
new_dat_of_clinidal %>% 
#dat_of_clinical %>% 
  select(paper_BRCA_Subtype_PAM50,subtype_PAM50_classifier,
         patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc,barcode) %>% 
  mutate(ImHi=map2_chr(er_status_by_ihc,pr_status_by_ihc,~(paste(.x,.y,sep = ", ")))) %>% 
  mutate(ImHi = map2_chr(ImHi,her2_status_by_ihc,~(paste(.x,.y,sep = ", ")))) %>%  group_by(patient) %>% 
  nest() %>%  
  mutate(number = map_dbl(data,nrow)) %>% filter(number>1) %>% 
  mutate(subtype = map_chr(data,~(paste0(.$cluster_pearson,collapse = ", ")))) %>% 
  mutate(ImHi = map(data,"ImHi"))

#exclude dupliated samples
#to explore relationship between clustering results and immunohistology subtype
dat_for_imhi_subtype <- new_dat_of_clinidal %>% filter(!duplicated(patient)) %>% 
  select(paper_BRCA_Subtype_PAM50,subtype_PAM50_classifier,
         patient,cluster_euclidean,cluster_pearson_centering,cluster_pearson,
         er_status_by_ihc,pr_status_by_ihc,her2_status_by_ihc,barcode) %>% 
  mutate(ImHi=map2_chr(er_status_by_ihc,pr_status_by_ihc,~(paste(.x,.y,sep = ", ")))) %>% 
  mutate(ImHi = map2_chr(ImHi,her2_status_by_ihc,~(paste(.x,.y,sep = ", "))))

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

dat_with_imhi_subtype <- dat_for_imhi_subtype %>% 
  mutate(imhi_subtype = map_chr(ImHi,eval_imhi_subtype,imhi_subtype_def))


# New prototype is better -------------------------------------------------

#check ER/PR/Her status 
#in euclidean clustering
dat_with_imhi_subtype %>% 
  filter(!grepl("Not|Indeterminate|Equivocal",ImHi)) %>% 
  group_by(ImHi) %>% count(cluster_euclidean) %>% group_by(cluster_euclidean) %>% 
  summarise(a=n/sum(n),b=ImHi,n=n)


#in PAM50 subtyping res
dat_with_imhi_subtype %>% 
  filter(!grepl("Not|Indeterminate|Equivocal",ImHi)) %>% 
  group_by(ImHi) %>% count(subtype_PAM50_classifier) %>% group_by(subtype_PAM50_classifier) %>% 
  summarise(a=n/sum(n),b=ImHi,n=n)

#in pearson_centering clustering
dat_with_imhi_subtype %>% 
  filter(!grepl("Not|Indeterminate|Equivocal",ImHi)) %>% 
  group_by(ImHi) %>% count(cluster_pearson_centering) %>% group_by(cluster_pearson_centering) %>% 
  summarise(a=n/sum(n),b=ImHi,n=n)

#combined_analysis with mutation_info
load("20221215_BRCA_mutation.RData")
mutation_info <- maf %>% 
  filter(!One_Consequence%in%"synonymous_variant") %>% group_by(Tumor_Sample_Barcode,Hugo_Symbol) %>% 
  summarise(mutated_number=n()) %>% pivot_wider(names_from = Hugo_Symbol,values_from = mutated_number) %>% 
  mutate(patient=str_replace_all(Tumor_Sample_Barcode,"(TCGA)-(..)-(....).*","\\1-\\2-\\3")) %>% 
  relocate(patient,.after = Tumor_Sample_Barcode)

dat_with_imhi_sub_mut_num <- dat_with_imhi_subtype %>% 
  left_join(mutation_info,by = "patient")

#transform NA to 0, because NA represents no mutation and compare
dat_with_imhi_sub_mut <- dat_with_imhi_sub_mut_num %>% select(1:14) %>% 
  bind_cols(dat_with_imhi_sub_mut_num %>% select(-1:-14) %>% map_df(~(ifelse(is.na(.),0,1))))


#oncoplot by maftools
maftools_maf <- read.maf(maf = maf %>% 
                           mutate(Tumor_Sample_Barcode=str_replace_all(Tumor_Sample_Barcode,
                                                                       "(TCGA)-(..)-(....).*","\\1-\\2-\\3")) %>% 
                           filter(Tumor_Sample_Barcode%in%dat_with_imhi_subtype$patient),
                         clinicalData = dat_with_imhi_subtype %>% rename(Tumor_Sample_Barcode=patient))
oncoplot(
  maf = maftools_maf, #genes = c("BRCA1","BRCA2"),
  clinicalFeatures = c('cluster_euclidean','cluster_pearson_centering','cluster_pearson','subtype_PAM50_classifier'),
  top = 10,
  sortByAnnotation = TRUE
)

#seems the results are similar between cluster_euclidean and cluster_pearson_centering and subtype_PAM50_classifier



# If the clustering results are the same for array and RNA-seq data --------

load("20221216_BRCA_microarray_expression.RData")
assays(BRCA_microarray)











