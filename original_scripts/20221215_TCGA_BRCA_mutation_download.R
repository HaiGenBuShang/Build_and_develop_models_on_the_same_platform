library(TCGAbiolinks)
library(SummarizedExperiment)
library(reshape2)
library(data.table)
library(tidyverse)

setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

#BRCA
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation"
)
GDCdownload(query)
maf <- GDCprepare(query)
save(maf,file = "20221215_BRCA_mutation.RData")



# load("20221215_BRCA_mutation.RData")
# #mutation to R
# maf %>% group_by(Tumor_Sample_Barcode) %>% select(Tumor_Sample_Barcode,Chromosome,
#                                                   Start_Position,End_Position,
#                                                   Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,
#                                                   Variant_Type,Variant_Classification)
# mutation_dat <- maf %>%
#   select(Tumor_Sample_Barcode,Chromosome,
#          Start_Position,End_Position,
#          Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,
#          Variant_Type,Variant_Classification) %>%
#   group_by(Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position) %>%
#   summarise(mutation=paste(Chromosome,Start_Position,End_Position,
#                            Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2,
#                            sep = "_",collapse = "")) %>% data.table() %>%
#   dcast(Tumor_Sample_Barcode~mutation,fun.aggregate = length)



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

save(BRCA_array_expr,BRCA_array_clin,file = "20221216_BRCA_microarray_expression.RData")



