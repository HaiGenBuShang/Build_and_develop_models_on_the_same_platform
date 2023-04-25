#****************************************************************************
#*Supplementary Figure 4 in manuscript
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
library(biomaRt)
#library(psych)
library(seqc)
library(dendextend)



#*******************************
#*SEQC clustering with TCGA BRCA
#*******************************
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


tumor_expr <- apply(tumor_expr,2,function(x){
  (x-mean(x))/sd(x)
})


SEQC_ILM_BGI_FPKM_entrez <- apply(ILM_refseq_gene_BGI[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_BGI[,3])

SEQC_ILM_BGI_FPKM_symbol <- apply(SEQC_ILM_BGI_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_BGI[,2])

SEQC_expr_BGI <- log2(SEQC_ILM_BGI_FPKM_symbol+1)


SEQC_ILM_CNL_FPKM_entrez <- apply(ILM_refseq_gene_CNL[,-1:-4],2,function(x,y){
  x/((y/1000)*(sum(x)/10e5))
},y=ILM_refseq_gene_CNL[,3])

SEQC_ILM_CNL_FPKM_symbol <- apply(SEQC_ILM_CNL_FPKM_entrez,2,function(x,y){
  tapply(x,y,mean)
},y=ILM_refseq_gene_CNL[,2])

SEQC_expr_CNL <- log2(SEQC_ILM_CNL_FPKM_symbol+1)





dat_1 <- SEQC_expr_CNL[rownames(SEQC_expr_CNL),]
colnames(dat_1) <- paste0(colnames(dat_1),"_CNL")
dat_2 <- SEQC_expr_BGI[rownames(SEQC_expr_BGI),]
colnames(dat_2) <- paste0(colnames(dat_2),"_BGI")

TCGA_and_SEQC_data <- cbind(dat_1,dat_2)
sample_color <- rep(c("#377EB8","#4DAF4A"),
                    times=c(ncol(dat_1),ncol(dat_2)))
names(sample_color) <- colnames(TCGA_and_SEQC_data)

clustering_res <- hclust(d = dist(t(TCGA_and_SEQC_data)),method = "ward.D") %>% as.dendrogram()

pdf("../Results/Figure/Supplementary_Figure_4.pdf",width = 20,height = 10)
clustering_res %>%  
  set("labels_cex",0.1) %>% 
  assign_values_to_leaves_edgePar(value = ifelse(grepl("TCGA",labels(clustering_res)),"#E41A1C","#377EB8"),
                                  edgePar = "col") %>%
  set("labels_colors",value = sample_color[labels(clustering_res)]) %>% 
  plot(main="Original Expression data",cex.main=2)

legend("topright",legend = c("SEQC BGI","SEQC CNL"),
       col=c("#377EB8","#4DAF4A"),
       lwd=1,bty = "n",cex = 2)
dev.off()
