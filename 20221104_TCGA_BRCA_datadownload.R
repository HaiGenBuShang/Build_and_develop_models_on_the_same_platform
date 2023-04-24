library(TCGAbiolinks)
library(SummarizedExperiment)

set.seed(12345678)
#down load TCGA BRCA expressions
setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")

# ##try TCGA-BRCA lagacy data 
# microarray data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "AgilentG4502A_07_3",
  experimental.strategy = "Gene expression array",
  legacy = TRUE
)

GDCdownload(query)
expr <- GDCprepare(query)

save(expr,file = "TCGA_BRCA_array_expr_GDCprepare.RData")