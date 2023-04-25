#****************************************************************************
#*all figures in manuscript by Fig. 1 and supplementary Figure 2
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


# Figure 1 ----------------------------------------------------------------

rm(list = ls())
# RNA-seq data
set.seed(12345678)
#down load TCGA BRCA expressions
setwd("./Data/")

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


#for anno data
wd <- getwd()
setwd("./")
TCGA_symbol_and_entrez <- read.table("TCGA_gene_symbol_and_entrez_ID.txt",
                                     header = TRUE,sep = "\t",stringsAsFactors = FALSE)
setwd(wd)

#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
#scale RNA-seq data
tumor_expr <- t(apply(tumor_expr,1,scale))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)


tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


#array data 
tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
tumor_expr <- t(apply(tumor_expr,1,scale))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                             annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                stringsAsFactors = FALSE),
                                                             do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_array_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                      annot = anno,
                                                      do.mapping = FALSE)

clustering_dat_for_plot <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype,
                                 PAM50_subtyping_res_RNA_seq_AIMS$subtype[,1],
                                 PAM50_subtyping_res_array_AIMS$subtype[,1]) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2,RNA_seq_AIMS=3,Array_AIMS=4)


pdf("../Results/Figure/Figure_1b.pdf",width = 10,height = 5)
clustering_dat_for_plot %>% arrange(RNA_seq_PAM50,Array_PAM50,RNA_seq_AIMS,Array_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_PAM50","Array_PAM50","RNA_seq_AIMS","Array_AIMS")[4:1])) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank(),axis.text.y = element_text(size=15),legend.position = "bottom") +
  labs(x='',y='',fill='Subtypes') 
dev.off()





##*********************
##*median centered data
##*********************
#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
#scale RNA-seq data
tumor_expr <- t(apply(tumor_expr,1,function(x)x-median(x)))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)


tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


#array data 
tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
tumor_expr <- t(apply(tumor_expr,1,function(x)x-median(x)))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                             annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                stringsAsFactors = FALSE),
                                                             do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_array_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                      annot = anno,
                                                      do.mapping = FALSE)

clustering_dat_for_plot <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype,
                                 PAM50_subtyping_res_RNA_seq_AIMS$subtype[,1],
                                 PAM50_subtyping_res_array_AIMS$subtype[,1]) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2,RNA_seq_AIMS=3,Array_AIMS=4)


pdf("../Results/Figure/Supplementary_Figure_2b_median_centered.pdf",width = 10,height = 5)
clustering_dat_for_plot %>% arrange(RNA_seq_PAM50,Array_PAM50,RNA_seq_AIMS,Array_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_PAM50","Array_PAM50","RNA_seq_AIMS","Array_AIMS")[4:1])) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank(),axis.text.y = element_text(size=15),legend.position = "bottom") +
  labs(x='',y='',fill='Subtypes') 
dev.off()

##*********************
##*mean centered data
##*********************
#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
#scale RNA-seq data
tumor_expr <- t(apply(tumor_expr,1,function(x)x-mean(x)))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)


tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


#array data 
tumor_expr <- array_expr_for_clust[,intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
tumor_expr <- t(apply(tumor_expr,1,function(x)x-mean(x)))
colnames(tumor_expr) <- intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))

data("pam50")
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                             annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                stringsAsFactors = FALSE),
                                                             do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_array_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                      annot = anno,
                                                      do.mapping = FALSE)

clustering_dat_for_plot <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype,
                                 PAM50_subtyping_res_RNA_seq_AIMS$subtype[,1],
                                 PAM50_subtyping_res_array_AIMS$subtype[,1]) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2,RNA_seq_AIMS=3,Array_AIMS=4)


pdf("../Results/Figure/Supplementary_Figure_2a_median_centered.pdf",width = 10,height = 5)
clustering_dat_for_plot %>% arrange(RNA_seq_PAM50,Array_PAM50,RNA_seq_AIMS,Array_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_PAM50","Array_PAM50","RNA_seq_AIMS","Array_AIMS")[4:1])) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank(),axis.text.y = element_text(size=15),legend.position = "bottom") +
  labs(x='',y='',fill='Subtypes') 
dev.off()


#Original value subtyping using PAM50 and AIMS
RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]
array_expr_PAM50 <- array_expr_for_clust[,
                                         intersect(colnames(RNA_seq_expr_for_clust),colnames(array_expr_for_clust))]

#RNA_seq data
tumor_expr <- RNA_seq_expr_PAM50
data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


#array data
tumor_expr <- array_expr_PAM50
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                             annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                stringsAsFactors = FALSE),
                                                             do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_array_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                      annot = anno,
                                                      do.mapping = FALSE)

clustering_dat_for_plot_2 <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype,
                                   PAM50_subtyping_res_RNA_seq_AIMS$subtype[,1],
                                   PAM50_subtyping_res_array_AIMS$subtype[,1]) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2,RNA_seq_AIMS=3,Array_AIMS=4)

pdf("../Results/Figure/Figure_1a.pdf",width = 10,height = 5)
clustering_dat_for_plot_2 %>% arrange(RNA_seq_PAM50,Array_PAM50,RNA_seq_AIMS,Array_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_PAM50","Array_PAM50","RNA_seq_AIMS","Array_AIMS")[4:1])) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank(),axis.text.y = element_text(size=15),legend.position = "bottom") +
  labs(x='',y='',fill='Subtypes')
dev.off()



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


clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D")

samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1)


samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n!=2)



##**************
##*New Figure 1c
##**************
consensus_samples_array_and_RNA_seq <- samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

(samples_of_5_classes <- cutree(clustering_res,k=5) %>% as.matrix() %>% as.data.frame() %>% 
    rownames_to_column("barcode") %>% rename(cluster_num=2) %>% mutate(cluster_num=as.character(cluster_num)) %>% 
    mutate(samples=str_replace_all(barcode,"_.*","")))


consensus_samples_of_5_classes <- samples_of_5_classes %>% 
  inner_join(consensus_samples_array_and_RNA_seq,by="samples") %>%
  rename(Subtype=2) %>% select(1,2) %>%   
  mutate(Clustering="Clustered by samples") %>%
  full_join(samples_of_5_classes,by="barcode") %>% 
  mutate(Clustering=ifelse(is.na(Clustering),"Not clustered by samples","Clustered by samples")) #%>% 
#mutate(Subtype=ifelse(is.na(Subtype),"Not Available",Subtype))



annotation_dat <- consensus_samples_of_5_classes %>% select(1,2,3) %>% rename(Subtype=2) %>% 
  rename(Subtype=2,barcode=1)

annotation_col <- data.frame(annotation_dat,row.names = 1)
anno_color <- list(Subtype=set_names(RColorBrewer::brewer.pal(5,"Set1"),as.character(1:5)),
                   Clustering=set_names(c(RColorBrewer::brewer.pal(5,"Set1")[2:1]),
                                        c("Clustered by samples","Not clustered by samples")))

pdf("../Results/Figure/Figure_1c.pdf",width = 10,height = 5)
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),".*",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "ward.D",annotation_colors = anno_color,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))
dev.off()




#Original Kappa
caret::confusionMatrix(data = factor(clustering_dat_for_plot_2$RNA_seq_PAM50,
                                     levels = c("Basal","Her2","LumA","LumB","Normal")),
                       reference = factor(clustering_dat_for_plot_2$Array_PAM50,
                                          levels = c("Basal","Her2","LumA","LumB","Normal")))
caret::confusionMatrix(data = factor(clustering_dat_for_plot_2$RNA_seq_AIMS,
                                     levels = c("Basal","Her2","LumA","LumB","Normal")),
                       reference = factor(clustering_dat_for_plot_2$Array_AIMS,
                                          levels = c("Basal","Her2","LumA","LumB","Normal")))

#Sclaed Kappa
caret::confusionMatrix(data = factor(clustering_dat_for_plot$RNA_seq_PAM50,
                                     levels = c("Basal","Her2","LumA","LumB","Normal")),
                       reference = factor(clustering_dat_for_plot$Array_PAM50,
                                          levels = c("Basal","Her2","LumA","LumB","Normal")))
caret::confusionMatrix(data = factor(clustering_dat_for_plot$RNA_seq_AIMS,
                                     levels = c("Basal","Her2","LumA","LumB","Normal")),
                       reference = factor(clustering_dat_for_plot$Array_AIMS,
                                          levels = c("Basal","Her2","LumA","LumB","Normal")))


##select Basal samples to predicting
#RNA_seq data
tumor_expr <- RNA_seq_expr_for_clust[,(clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames()]
#scale RNA-seq data
tumor_expr <- t(apply(tumor_expr,1,scale))
colnames(tumor_expr) <- (clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames()

data("pam50")
PAM50_subtyping_res_RNA_seq_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                               annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                  stringsAsFactors = FALSE),
                                                               do.mapping = FALSE)


tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_RNA_seq_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                        annot = anno,
                                                        do.mapping = FALSE)


#array data 
tumor_expr <- array_expr_for_clust[,
                                   (clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames()]
tumor_expr <- t(apply(tumor_expr,1,scale))
colnames(tumor_expr) <- (clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames()

data("pam50")
PAM50_subtyping_res_array_PAM50 <- intrinsic.cluster.predict(sbt.model = pam50,data = t(tumor_expr),
                                                             annot = data.frame(Gene.Symbol=rownames(tumor_expr),
                                                                                stringsAsFactors = FALSE),
                                                             do.mapping = FALSE)
tumor_expr <- 2^tumor_expr
tumor_expr <- tumor_expr[!is.na(match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol)),]
anno <- TCGA_symbol_and_entrez[match(rownames(tumor_expr),TCGA_symbol_and_entrez$Gene.Symbol),]
PAM50_subtyping_res_array_AIMS <- molecular.subtyping(sbt.model = "AIMS",data = t(tumor_expr),
                                                      annot = anno,
                                                      do.mapping = FALSE)

clustering_dat_for_plot_3 <- cbind(PAM50_subtyping_res_RNA_seq_PAM50$subtype,PAM50_subtyping_res_array_PAM50$subtype,
                                   PAM50_subtyping_res_RNA_seq_AIMS$subtype[,1],
                                   PAM50_subtyping_res_array_AIMS$subtype[,1]) %>% 
  as.data.frame() %>% mutate(barcode=rownames(.)) %>% 
  rename(RNA_seq_PAM50=1,Array_PAM50=2,RNA_seq_AIMS=3,Array_AIMS=4)


pdf("../Results/Figure/Figure_1e.pdf",width = 10,height = 5)
clustering_dat_for_plot_3 %>% arrange(RNA_seq_PAM50,Array_PAM50,RNA_seq_AIMS,Array_AIMS) %>% 
  mutate(barcode=factor(barcode,levels = unique(barcode))) %>% 
  pivot_longer(cols = -barcode,names_to = "Clustering",values_to = "value") %>% 
  mutate(Clustering=factor(Clustering,levels = c("RNA_seq_PAM50","Array_PAM50","RNA_seq_AIMS","Array_AIMS")[4:1])) %>% 
  ggplot() +
  geom_tile(aes(x=barcode,y=Clustering,fill=factor(value))) +
  #scale_fill_viridis(discrete=FALSE) +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.background = element_blank(),axis.text.y = element_text(size=15),legend.position = "bottom") +
  labs(x='',y='',fill='Subtypes') 
dev.off()


#Basal sample kappa
caret::confusionMatrix(data = factor(clustering_dat_for_plot_3$RNA_seq_PAM50,
                                     levels = c("Basal","Her2","LumA","LumB","Normal")),
                       reference = factor(clustering_dat_for_plot_3$Array_PAM50,
                                          levels = c("Basal","Her2","LumA","LumB","Normal")))
caret::confusionMatrix(data = factor(clustering_dat_for_plot_3$RNA_seq_AIMS,
                                     levels = c("Basal","Her2","LumA","LumB","Normal")),
                       reference = factor(clustering_dat_for_plot_3$Array_AIMS,
                                          levels = c("Basal","Her2","LumA","LumB","Normal")))




# Basal sample clustering --------------------------------

RNA_seq_expr_PAM50 <- RNA_seq_expr_for_clust[,
                                             (clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames()]
array_expr_PAM50 <- array_expr_for_clust[,
                                         (clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames()]
intersect_genes <- intersect(rownames(RNA_seq_expr_PAM50),rownames(array_expr_PAM50))
RNA_seq_expr_PAM50 <- t(apply(RNA_seq_expr_PAM50[intersect_genes,],1,scale))
array_expr_PAM50 <- t(apply(array_expr_PAM50[intersect_genes,],1,scale))

combined_dat <- cbind(RNA_seq_expr_PAM50,array_expr_PAM50)
colnames(combined_dat) <- c(paste0((clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames(),
                                   "_RNAseq"),
                            paste0((clustering_dat_for_plot %>% filter(RNA_seq_PAM50=="Basal")) %>% rownames(),
                                   "_array"))



clustering_res <- hclust(d = dist(t(combined_dat)),method = "ward.D")
samples_of_classes <- cutree(clustering_res,k = ncol(RNA_seq_expr_PAM50)) %>% as.matrix() %>% as.data.frame() %>%
  mutate(barcode=rownames(.)) %>% rename(cluster_num=1) %>% mutate(cluster_num=as.character(cluster_num))


consensus_samples_array_and_RNA_seq <- samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

consensus_samples_array_and_RNA_seq <- samples_of_classes %>% mutate(samples=str_replace_all(barcode,"_.*","")) %>% 
  group_by(samples,cluster_num) %>% count() %>% filter(n==2)

(samples_of_5_classes <- cutree(clustering_res,k=5) %>% as.matrix() %>% as.data.frame() %>% 
    rownames_to_column("barcode") %>% rename(cluster_num=2) %>% mutate(cluster_num=as.character(cluster_num)) %>% 
    mutate(samples=str_replace_all(barcode,"_.*","")))


consensus_samples_of_5_classes <- samples_of_5_classes %>% 
  inner_join(consensus_samples_array_and_RNA_seq,by="samples") %>%
  rename(Subtype=2) %>% select(1,2) %>%   
  mutate(Clustering="Clustered by samples") %>%
  full_join(samples_of_5_classes,by="barcode") %>% 
  mutate(Clustering=ifelse(is.na(Clustering),"Not clustered by samples","Clustered by samples")) #%>% 


annotation_dat <- consensus_samples_of_5_classes %>% select(1,2,3) %>% rename(Subtype=2) %>% 
  rename(Subtype=2,barcode=1)

annotation_col <- data.frame(annotation_dat,row.names = 1) %>% select(-1)
anno_color <- list(Subtype=set_names(RColorBrewer::brewer.pal(5,"Set1"),as.character(1:5)),
                   Clustering=set_names(c(RColorBrewer::brewer.pal(5,"Set1")[2:1]),
                                        c("Clustered by samples","Not clustered by samples")))

pdf("../Results/Figure/Figure_1d.pdf",width = 10,height = 5)
pheatmap(mat = combined_dat,labels_row = rep('',nrow(combined_dat)),
         labels_col = str_replace_all(colnames(combined_dat),".*",""),
         annotation_col = annotation_col,fontsize_col = 2,clustering_method = "ward.D",annotation_colors = anno_color,
         color = colorRampPalette(c("green", "black", "firebrick3"))(101),breaks = seq(-2,2,length.out=101))
dev.off()







