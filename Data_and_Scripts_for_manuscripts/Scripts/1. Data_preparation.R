library(TCGAbiolinks)
library(SummarizedExperiment)
library(reshape2)
library(data.table)
library(tidyverse)

setwd("./Data/")

#BRCA
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation"
)
GDCdownload(query)
maf <- GDCprepare(query)
save(maf,file = "20221215_BRCA_mutation.RData")

# Prepare microarray data -------------------------------------------------
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "AgilentG4502A_07_3",
  legacy = TRUE
)
GDCdownload(query)

## The command below would produce wrong data,
## so the data need to be prepared by hand

BRCA_microarray <- GDCprepare(query)

#all microarray data
all_array_files <- list.files(path = "./GDCdata/TCGA-BRCA/legacy/Gene_expression/Gene_expression_quantification/",
                              pattern = "level3.data.txt",recursive = TRUE,full.names = TRUE)

array_expr <- lapply(all_array_files,function(x){
  read.table(x,stringsAsFactors = FALSE,sep = "\t",header = TRUE,quote = "",check.names = FALSE)[-1,]
})

BRCA_array_expr <- array_expr %>% reduce(left_join,by="Hybridization REF")
BRCA_array_expr <- BRCA_array_expr %>% data.frame(row.names = 1, check.names = FALSE)
BRCA_array_clin <- colData(BRCA_microarray)

# BRCA_array_expr_sample_order <- read.table("BRCA_array_expr_samples.txt",
#                                            header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
# BRCA_array_clin_sample_order <- read.table("BRCA_array_clin_samples.txt",
#                                            header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
# 
# BRCA_array_expr <- BRCA_array_expr[,BRCA_array_expr_sample_order]
# BRCA_array_clin <- BRCA_array_clin[BRCA_array_clin_sample_order,]

save(BRCA_array_expr,BRCA_array_clin,file = "20221216_BRCA_microarray_expression.RData")





# Prepare RNA-seq data -------------------------------------------------
#download TCGA expression data
#find TCGA project
TCGA_projects <- "TCGA-BRCA"

#query TCGA data by project
TCGA_query <- lapply(TCGA_projects,function(x){
  GDCquery(project = x,
           data.category = "Transcriptome Profiling",
           data.type = "Gene Expression Quantification",
           workflow.type = "STAR - Counts")
})
names(TCGA_query) <- "TCGA-BRCA"

#Download TCGA data by project
lapply(TCGA_query,GDCdownload)

#BRCA prepare
dir.create(path = "R_data/BRCA",recursive = TRUE)
expr <- GDCprepare(TCGA_query$`TCGA-BRCA`)
X1 <- assay(expr,"fpkm_unstrand")
A <- colData(expr) %>% as.data.frame()
gene_info <- rowData(expr) %>% as.data.frame()
X1_genename <- apply(X1,2,function(x){
  tapply(x,gene_info$gene_name,mean)
})

# X1_sample_order <- read.table("X1_samples.txt",
#                               header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
# A_sample_order <- read.table("A_samples.txt",
#                              header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
# 
# X1 <- X1[,X1_sample_order]
# X1_genename <- X1_genename[,X1_sample_order]
# A <- A[A_sample_order,]


save(X1,A,X1_genename,file = "R_data/BRCA/Expr_and_pheno.RData")



