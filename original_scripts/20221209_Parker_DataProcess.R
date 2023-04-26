#Parker data subtype
library(tidyverse)
setwd("/mnt/Miscrosoft/Shi_lab/Breast_cancer/TNBC/PAM50/Parker_data_microarray/txt_data/")



dt_anno <- read.table("GPL1390_annotation.txt",header = TRUE,stringsAsFactors = FALSE,fill = TRUE,sep = "\t",
                      skip = 1)
annot_info <- dt_anno[,c(1,7)]
write.table(annot_info,file = "GPL1390_annotation_simplified_2.txt",sep = "\t",col.names = TRUE,row.names = FALSE)



##parker_data
read_parker <- function(data_file,skip=83){
  dt <- read.table(data_file,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = skip,
                   fill = TRUE)
  dt <- dt[-nrow(dt),]
}
##map parker
map_parker <- function(expr_dt,annot_file,skip_anno=1){
  dt_anno <- read.table(annot_file,header = TRUE,stringsAsFactors = FALSE,fill = TRUE,sep = "\t",
                        skip = skip_anno,quote = "")
  expr_dt[,1] <- dt_anno[,2]
  expr_dt <- subset(expr_dt,ID_REF!="")
  expr_dt_mean_probe_exprs_for_each_gene <- apply(expr_dt[,-1],2,function(x){
    tapply(x,expr_dt[,1],mean)
  })
}

Parker_expr_GPL885 <- read_parker(
  data_file = "GSE10886-GPL885_series_matrix.txt",
  skip = 83)
Parker_expr_GPL885 <- map_parker(Parker_expr_GPL885,
                                 annot_file = paste0("GPL885_annotation_simplified.txt"))

Parker_expr_GPL887 <- read_parker(
  data_file = "GSE10886-GPL887_series_matrix.txt",
  skip = 85)
Parker_expr_GPL887 <- map_parker(Parker_expr_GPL887,
                                 annot_file = paste0("GPL887_annotation_simplified.txt"))


Parker_expr_GPL1390 <- read_parker(
  data_file = "GSE10886-GPL1390_series_matrix.txt",
  skip = 88)
Parker_expr_GPL1390 <- map_parker(Parker_expr_GPL1390,
                                 annot_file = paste0("GPL1390_annotation_simplified_2.txt"),skip_anno = 0)

save.image(file = "20221209_Parker_GPL_expr_data.RData")


Parker_clinical <- read.table("GSE10886-GPL1390_series_matrix.txt",header = TRUE,stringsAsFactors = FALSE,
                              sep = "\t",skip = 31,nrows = 55,check.names = FALSE)

Parker_clinical_GPL885 <- read.table("GSE10886-GPL885_series_matrix.txt",header = TRUE,stringsAsFactors = FALSE,
                                     sep = "\t",skip = 31,nrows = 50,check.names = FALSE)

Parker_clinical_GPL887 <- read.table("GSE10886-GPL887_series_matrix.txt",header = TRUE,stringsAsFactors = FALSE,
                                     sep = "\t",skip = 31,nrows = 52,check.names = FALSE)


Parker_clinical_GPL887 %>% slice(16) %>% pivot_longer(-1)
Parker_clinical %>% slice(16) %>% pivot_longer(-1)
Parker_clinical_GPL885 %>% slice(16) %>% pivot_longer(-1)

obtain_T_and_N <- function(clinical_data){
  sample_type <- clinical_data %>% slice(16) %>% pivot_longer(-1)
  sample_type_T_N <- sample_type %>% 
    filter(grepl("breast tumor",value,ignore.case = TRUE)|grepl("normal",value,ignore.case = TRUE))
  
  sample_type_T_N_2 <- sample_type %>% 
    filter(!(grepl("breast tumor",value,ignore.case = TRUE)|grepl("normal",value,ignore.case = TRUE))) %>% 
    mutate(sample_detail = map(name,~(select(clinical_data %>% filter(grepl("ch2",`!Sample_title`)),.)))) %>%
    mutate(Tumor = map_lgl(sample_detail,~(any(grepl("breast tumor",.,ignore.case = TRUE))))) %>% 
    mutate(Normal=map_lgl(sample_detail,~(any(grepl("normal",.,ignore.case = TRUE))))) %>% 
    mutate(losted=map2_lgl(Tumor,Normal,~(any(c(.x,.y))))) %>% 
    mutate(value = as.character(ifelse(Tumor,"Breast tumor",ifelse(Normal,"Normal",NA)))) #%>% select(1:3)
  bind_rows(sample_type_T_N,sample_type_T_N_2)
  
  #sample_type_T_N_2 %>% filter(!losted)
}


Parker_clinical_GPL1390_T_N <- obtain_T_and_N(clinical_data = Parker_clinical)

Parker_clinical_GPL1390_T_N <- Parker_clinical_GPL1390_T_N %>% 
  mutate(value = ifelse(name=="BR01-0125-LumB","Breast tumor",value)) %>% 
  mutate(value = ifelse(name=="UB-120","Breast tumor",value)) %>% 
  mutate(value = ifelse(name=="UB131","Breast tumor",value)) %>% filter(!is.na(value))


Parker_clinical_GPL885_T_N <- obtain_T_and_N(clinical_data = Parker_clinical_GPL885)
Parker_clinical_GPL887_T_N <- obtain_T_and_N(clinical_data = Parker_clinical_GPL887)

bind_rows(Parker_clinical_GPL1390_T_N,Parker_clinical_GPL885_T_N,Parker_clinical_GPL887_T_N) %>% count(value)


Parker_clinical_GPL1390_T_N %>% mutate(list=map_lgl(sample_detail,is.null)) %>% filter(!list)














