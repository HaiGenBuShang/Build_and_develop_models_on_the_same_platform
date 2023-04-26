library(breastCancerMAINZ)
library(sigclust)
library(tidyverse)
library(sigclust2)

setwd("/mnt/Miscrosoft/Brease_Cancer_subtyping/Data/")


data("mainz")
g_expr <- exprs(mainz)

#probe expression to gene symbol expression by mean
g_anno <- fData(mainz) %>% select(probe,Gene.symbol)

g_expr_symbol <- apply(g_expr,2,function(x,y){
  tapply(x,INDEX = y,mean,na.rm=TRUE)
},y=g_anno$Gene.symbol)

#intrinsic gene expression
intrinsic_genes <- read.table("2009_JCO_intrinsic_genes_S_table_5.txt",header = FALSE,sep = "\t",
                              stringsAsFactors = FALSE)

intrinsic_g_expr <- g_expr_symbol[intersect(rownames(g_expr_symbol),intrinsic_genes[,1]),]
intrinsic_expr_median_centered <- apply(intrinsic_g_expr,1,function(x)x-mean(x))


cluster_res <- shc(intrinsic_expr_median_centered,metric="cor", linkage="average",null_alg = "2means")











