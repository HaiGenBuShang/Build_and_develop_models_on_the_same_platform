library(genefu)
library(TCGA2STAT)
library(tidyverse)
##obtain PAM50 genes
ob_PAM50_genes <- function(expr_dt){
  PAM50_centroids <- read.table(
    "./Data/pam50_centroids.txt",
    sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  shared_genes <- intersect(rownames(expr_dt),PAM50_centroids[,1])
  expr_dt[shared_genes,]
}


message("Please note that if TCGA2STAT::getTCGA throw errors
you can using TCGAbiolinks downloaded data and the results would be a little difference
but not affect the conclusion
The TCGA2STAT might be throw many errors and it would consume much time to repair them")



set.seed(seed = 12345678)

TCGA_BRCA_expr <- TCGA2STAT::getTCGA(disease = "BRCA",clinical = TRUE)

#dim(TCGA_BRCA_expr$dat)
TCGA_BRCA_expr$dat <- log2(TCGA_BRCA_expr$dat+0.01)
TCGA_BRCA_expr$dat <- TCGA_BRCA_expr$dat[,as.numeric(substr(colnames(TCGA_BRCA_expr$dat),14,15))<=1]
colnames(TCGA_BRCA_expr$dat) <- gsub("-","_",substr(colnames(TCGA_BRCA_expr$dat),1,15))
TCGA_BRCA_expr$merged.dat[,1] <- gsub("-","_",TCGA_BRCA_expr$merged.dat[,1])
TCGA_BRCA_expr$merged.dat$OS <- TCGA_BRCA_expr$merged.dat$OS/365


TCGA_PAM50_gene_expr <- ob_PAM50_genes(TCGA_BRCA_expr$dat)

median_expr <- t(apply(TCGA_PAM50_gene_expr,1,function(x)x - median(x)))
data("pam50")

cor(median_expr[,5],pam50$centroids[rownames(median_expr),],method = "spearman")
cor(TCGA_PAM50_gene_expr[,5],pam50$centroids[rownames(TCGA_PAM50_gene_expr),],method = "spearman")



pdf("./Results/Figure/Supplementary_Figure_1.pdf",width = 6.5,height = 6.5)
par(mai=c(1.52,1.32,1.32,0.92))
plot(TCGA_PAM50_gene_expr[,5],pam50$centroids[rownames(TCGA_PAM50_gene_expr),5],
     xlim=c(min(c(TCGA_PAM50_gene_expr[,5],median_expr[,5])),max(c(TCGA_PAM50_gene_expr[,5],median_expr[,5]))),
     xlab="Gene expression",ylab="Centroids",main="Spearman Correlation Before and\nAfter Median-Centering Normalization",
     col="#377EB8",pch=16,font.lab=2,cex.lab=1.4,cex.main=1.4,cex.axis=1.3)
fit <- lm(pam50$centroids[rownames(TCGA_PAM50_gene_expr),5]~TCGA_PAM50_gene_expr[,5])
abline(fit,col="#377EB8",lwd=2)

points(median_expr[,5],pam50$centroids[rownames(median_expr),5],col="#E41A1C",pch=16)
fit <- lm(pam50$centroids[rownames(median_expr),5]~median_expr[,5])
abline(fit,col="#E41A1C",lwd=2)

legend("bottomright",bty = "n",legend = c("Before","After"),lwd = 2,col = c("#377EB8","#E41A1C"),cex=1.5)
dev.off()












#***********************************
#*Using TCGAbiolinks downloaded data
#***********************************

load("./Data/R_data/BRCA/Expr_and_pheno.RData")
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

TCGA_PAM50_gene_expr <- ob_PAM50_genes(tumor_expr)
median_expr <- t(apply(TCGA_PAM50_gene_expr,1,function(x)x - median(x)))

data("pam50")
cor(median_expr[,6],pam50$centroids[rownames(median_expr),],method = "spearman") %>% abs() %>% which.max()
cor(TCGA_PAM50_gene_expr[,6],pam50$centroids[rownames(TCGA_PAM50_gene_expr),],method = "spearman")%>% abs() %>% which.max()

pdf("./Results/Figure/Supplementary_Figure_1.pdf",width = 6.5,height = 6.5)
par(mai=c(1.52,1.32,1.32,0.92))
plot(TCGA_PAM50_gene_expr[,6],pam50$centroids[rownames(TCGA_PAM50_gene_expr),5],
     xlim=c(min(c(TCGA_PAM50_gene_expr[,6],median_expr[,6])),max(c(TCGA_PAM50_gene_expr[,6],median_expr[,6]))),
     xlab="Gene expression",ylab="Centroids",main="Spearman Correlation Before and\nAfter Median-Centering Normalization",
     col="#377EB8",pch=16,font.lab=2,cex.lab=1.4,cex.main=1.4,cex.axis=1.3)
fit <- lm(pam50$centroids[rownames(TCGA_PAM50_gene_expr),5]~TCGA_PAM50_gene_expr[,6])
abline(fit,col="#377EB8",lwd=2)

points(median_expr[,6],pam50$centroids[rownames(median_expr),5],col="#E41A1C",pch=16)
fit <- lm(pam50$centroids[rownames(median_expr),5]~median_expr[,6])
abline(fit,col="#E41A1C",lwd=2)

legend("bottomright",bty = "n",legend = c("Before","After"),lwd = 2,col = c("#377EB8","#E41A1C"),cex=1.5)
dev.off()


