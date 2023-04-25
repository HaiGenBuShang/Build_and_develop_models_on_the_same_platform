#**********************************************************************************************
#******Reproduce PAM results on web http://statweb.stanford.edu/~tibs/PAM/comparison.html******
#**********************************************************************************************
library(sda)
library(pamr)
library(multtest)
library(HiDimDA)

#Define functions
produce_train_test_set <- function(expr_with_label,train_proportion=2/3){
  cat("The expr_with_label object, a list, must have two components x and y,\nx contain expression data and y contain labels!\n")
  if(ncol(expr_with_label$x)!=length(expr_with_label$y))
    stop("the numer of columns in x must equal to number of element in y\nand each column correspond to each element")
  expr_with_label$y <- as.character(expr_with_label$y)
  sample_n <- table(expr_with_label$y)
  sample_n_in_train <- round(sample_n*train_proportion,digits = 0)
  index <- 1:length(expr_with_label$y)
  expr_with_label$y <- setNames(expr_with_label$y,index)
  splited_samples <- split(expr_with_label$y,f = expr_with_label$y)
  
  train_sample_index <- mapply(function(x,y){
    names(sample(x = x,size = y))
  },x = splited_samples,y = sample_n_in_train,SIMPLIFY = FALSE)
  
  test_sample_index <- mapply(function(x,y){
    setdiff(names(x),y)
  },x = splited_samples, y = train_sample_index,SIMPLIFY = FALSE)
  
  train_sample_index <- as.integer(unlist(train_sample_index))
  test_sample_index <- as.integer(unlist(test_sample_index))
  train_set <- list(x=expr_with_label$x[,train_sample_index],y=expr_with_label$y[train_sample_index])
  test_set <- list(x=expr_with_label$x[,test_sample_index],y=expr_with_label$y[test_sample_index])
  
  if(is.null(rownames(train_set$x))){
    rownames(train_set$x) <- paste0("gene_",1:nrow(train_set$x))
  }
  if(is.null(rownames(test_set$x))){
    rownames(test_set$x) <- paste0("gene_",1:nrow(test_set$x))
  }
  
  list(train_set=train_set,test_set=test_set)
}


# PAM_test_error <- function(train_and_test_dataset,auto_pam_delt=FALSE,n_threshold=100){
#   #get train set and test set 
#   train_set <- train_and_test_dataset$train_set
#   test_set <- train_and_test_dataset$test_set
#   #using training set to train model
#   training_model <- pamr.train(train_set,n.threshold = n_threshold)
#   
#   cross_validation <- pamr.cv(fit = training_model,data = train_set,nfold = 10)
#   # 
#   # pamr.plotcv(cross_validation)
#   
#   #determine delt value of pam
#   #you can choose manully or determine it by choosing the mean of the lowest prediction error rate and C-V error rate
#   if(auto_pam_delt){
#     #obtain the delt which give the lowest training error
#     errors <- setNames(training_model$errors,1:length(training_model$errors))
#     errors <- errors[length(training_model$errors):1]
#     min_error_from_last <- as.integer(names(which.min(errors)))
#     min_error_delt <- training_model$threshold[min_error_from_last]
#     # pam_delt <- min_error_delt
#     
#     #obtain the delt which give the lowest Cross-validation error
#     CV_errors <- setNames(cross_validation$error,1:length(cross_validation$error))
#     CV_errors <- CV_errors[length(cross_validation$error):1]
#     min_CV_error_from_last <- as.integer(names(which.min(CV_errors)))
#     min_CV_error_delt <- cross_validation$threshold[min_CV_error_from_last]
#     # pam_delt <- min_CV_error_delt
#     # 
#     # 
#     # print(min_error_delt)
#     # print(min_CV_error_delt)
#     pam_delt <- mean(c(min_error_delt,min_CV_error_delt))
#   }else{
#     print(training_model)
#     print(cross_validation)
#     pam_delt <- as.numeric(readline(prompt = "which delt you choose? \n"))
#   }
#   test_predicted <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt)
#   test_predicted <- as.character(test_predicted)
#   # print(pam_delt)
#   sum(test_predicted!=test_set$y)/length(test_predicted)
# }


PAM_test_error <- function(train_and_test_dataset,auto_pam_delt=FALSE,n_threshold=100,
                           auto_delt_by=c("min_training_error","min_CV_error","mean_of_min_traing_and_CV_error")){
  #get train set and test set 
  train_set <- train_and_test_dataset$train_set
  test_set <- train_and_test_dataset$test_set
  #using training set to train model
  training_model <- pamr.train(train_set,n.threshold = n_threshold)
  
  cross_validation <- pamr.cv(fit = training_model,data = train_set,nfold = 10)
  ##cross_validation <- pamr.cv_with_Cor(fit = training_model,data = train_set,nfold = 10)
  # 
  # pamr.plotcv(cross_validation)
  
  #determine delt value of pam
  #you can choose manully or determine it by choosing the mean of the lowest prediction error rate and C-V error rate
  if(auto_pam_delt){
    #obtain the delt which give the lowest training error
    errors <- setNames(training_model$errors,1:length(training_model$errors))
    errors <- errors[length(training_model$errors):1]
    min_error_from_last <- as.integer(names(which.min(errors)))
    min_error_delt <- training_model$threshold[min_error_from_last]
    
    #obtain the delt which give the lowest training error
    CV_errors <- setNames(cross_validation$error,1:length(cross_validation$error))
    CV_errors <- CV_errors[length(cross_validation$error):1]
    min_CV_error_from_last <- as.integer(names(which.min(CV_errors)))
    min_CV_error_delt <- cross_validation$threshold[min_CV_error_from_last]
    
    switch(auto_delt_by,
           min_training_error={
             pam_delt <- min_error_delt
           },
           min_CV_error={
             pam_delt <- min_CV_error_delt
           },
           mean_of_min_traing_and_CV_error={
             pam_delt <- mean(c(min_error_delt,min_CV_error_delt))
           },
           stop("auto_delt_by arg must specified as one of (min_training_error min_CV_error mean_of_min_traing_and_CV_error)"))
  }else{
    print(training_model)
    print(cross_validation)
    pam_delt <- as.numeric(readline(prompt = "which delt you choose? \n"))
  }
  test_predicted <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt)
  test_predicted <- as.character(test_predicted)
  # print(pam_delt)
  sum(test_predicted!=test_set$y)/length(test_predicted)
}




n_train_and_test <- function(expr_with_label,n_times,train_proportion,auto_pam_delt,n_threshold,
                             ...){
                             #auto_delt_by=c("min_training_error","min_CV_error","mean_of_min_traing_and_CV_error"),...){
  testing_errors <- c()
  i <- 1
  while(i <= n_times){
    message(paste0("The ",i,"th repeat!"))
    train_and_test_dat <- produce_train_test_set(expr_with_label = expr_with_label,train_proportion = train_proportion)
    test_error <- PAM_test_error(train_and_test_dataset = train_and_test_dat,
                                 auto_pam_delt = auto_pam_delt,n_threshold = n_threshold,...)
    testing_errors <- c(testing_errors,test_error)
    i <- i+1
  }
  setNames(testing_errors,paste0("The ",1:n_times,"st errors"))
}


produce_dat_from_res_file <- function(res_file,max_d_min=5,max_m_min=500,min_signal=10,max_signal=16000){
  ##read data
  res_data <- read.table(res_file,header = FALSE,skip = 3,stringsAsFactors = FALSE,sep = "\t",quote = "",comment.char = "")
  headers <- read.table(res_file,header = FALSE,stringsAsFactors = FALSE,sep = "\t",nrows = 1)
  colnames(res_data) <- t(headers[1,])
  res_data <- res_data[,seq(3,ncol(res_data),by=2)]
  
  ##gene filtering
  res_data <- as.matrix(res_data)
  res_data <- res_data[,!apply(res_data,2,function(x)all(is.na(x)))]
  
  #assign min signal thres to low signal and max signal thres to high signal
  res_data[res_data<min_signal] <- min_signal
  res_data[res_data>max_signal] <- max_signal
  
  min_v <- apply(res_data,1,min,na.rm=TRUE)
  max_v <- apply(res_data,1,max,na.rm=TRUE)
  
  #two gene filter thres
  #one: max divide min > thres (e.g. 5)
  #two: max minus min > thres (eg. 500)
  filter_1 <- (max_v/min_v) > max_d_min
  filter_2 <- (max_v - min_v) > max_m_min
  
  res_data <- log10(res_data[filter_1&filter_2,])
  apply(res_data,2,function(x)(x-mean(x))/sd(x))
}
#res_file <- "d:/Shi_lab/Breast_cancer/New_PAM50_classifier/Data/data_summary/Dataset_A_multiple_tumor_samples.res"



#****************************************************************
#******Assessment PAM classification method and correlation******
#****************************************************************
library(sda)
library(pamr)
library(multtest)
library(HiDimDA)

#Define functions
# compare_pam_and_correlation <- function(train_and_test_dataset,auto_pam_delt=FALSE,n_threshold=100){
#   library(caret)
#   ##Pamr predict
#   #get train set and test set 
#   train_set <- train_and_test_dataset$train_set
#   test_set <- train_and_test_dataset$test_set
#   
#   if(is.null(rownames(train_set$x))){
#     rownames(train_set$x) <- paste0("gene_",1:nrow(train_set$x))
#   }
#   
#   if(is.null(rownames(test_set$x))){
#     rownames(test_set$x) <- paste0("gene_",1:nrow(test_set$x))
#   }
#   
#   #using training set to train model
#   training_model <- pamr.train(train_set,n.threshold = n_threshold)
#   
#   cross_validation <- pamr.cv(fit = training_model,data = train_set,nfold = 10)
#   # 
#   # pamr.plotcv(cross_validation)
#   
#   #determine delt value of pam
#   #you can choose manully or determine it by choosing the mean of the lowest prediction error rate and C-V error rate
#   if(auto_pam_delt){
#     #obtain the delt which give the lowest training error
#     # errors <- setNames(training_model$errors,1:length(training_model$errors))
#     # errors <- errors[length(training_model$errors):1]
#     # min_error_from_last <- as.integer(names(which.min(errors)))
#     # min_error_delt <- training_model$threshold[min_error_from_last]
#     # pam_delt <- min_error_delt
#     
#     #obtain the delt which give the lowest Cross-validation error
#     CV_errors <- setNames(cross_validation$error,1:length(cross_validation$error))
#     CV_errors <- CV_errors[length(cross_validation$error):1]
#     min_CV_error_from_last <- as.integer(names(which.min(CV_errors)))
#     min_CV_error_delt <- cross_validation$threshold[min_CV_error_from_last]
#     pam_delt <- min_CV_error_delt
#     # 
#     # 
#     # print(min_error_delt)
#     # print(min_CV_error_delt)
#     # pam_delt <- mean(c(min_error_delt,min_CV_error_delt))
#   }else{
#     print(training_model)
#     print(cross_validation)
#     pam_delt <- as.numeric(readline(prompt = "which delt you choose? \n"))
#   }
#   test_predicted <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt)
#   train_predicted <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt)
#   
#   ##correlation predicted
#   #assign geneid and genenames for training model,
#   #otherwise the pamr.listgenes would would throw error
#   train_set$geneid <- 1:nrow(train_set$x) #make the $geneid components in train_set the same as training_model$gene.subset
#   #if(!is.null(rownames(train_set$x))){
#     train_set$genenames <- rownames(train_set$x) #make the $genenames components in train_set the same as training_model$gene.subset
#   #}else{
#   #  train_set$genename <- as.character(1:nrow(train_set$x))
#   #}
#   
#   #pamr centroids
#   pamr_centroids <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "centroid")
#   pamr_survived_genes <- pamr.listgenes(fit = training_model,data = train_set,threshold = pam_delt,genenames = TRUE)[,"name"]
#   # pamr_survived_genes <- pamr.listgenes(fit = training_model,data = train_set,threshold = pam_delt,genenames = TRUE)[,"id"]
#   pamr_centroids <- pamr_centroids[pamr_survived_genes,]
#   
#   #correlation
#   pamr_centroids_test_cor <- cor(x = test_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
#   pamr_centroids_train_cor <- cor(x = train_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
#   
#   test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
#   train_cor_predicted <- apply(pamr_centroids_train_cor,1,function(x)names(which.max(x)))
#   
#   predict_res<-list(train_set=data.frame(pam_predict=train_predicted,
#                                          cor_predict=factor(train_cor_predicted,levels = levels(train_predicted)),
#                                          stringsAsFactors = FALSE),#set both column of data.frame as factor, so that
#                     #there would be the same columns in the table() function output, 
#                     #otherwise there would be an error in caret::confusionMatrix
#                     #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
#                     #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
#                     #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
#                     test_set=data.frame(pam_predict=test_predicted,
#                                         cor_predict=factor(test_cor_predicted,levels = levels(test_predicted)),
#                                         stringsAsFactors = FALSE)#set both column of data.frame as factor, so that
#                     #there would be the same columns in the table() function output, 
#                     #otherwise there would be an error in caret::confusionMatrix
#                     #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
#                     #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
#                     #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
#   )
#   kappas <- sapply(predict_res,function(x){
#     confusionMatrix(table(x))$overall["Kappa"]
#   })
#   
#   list(kappas=kappas,predict_res=predict_res,train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor,
#        centroid=pamr_centroids,delt=pam_delt)
# }


# compare_pam_and_correlation <- function(train_and_test_dataset,auto_pam_delt=FALSE,n_threshold=100,
#                                         auto_delt_by=c("min_training_error","min_CV_error","mean_of_min_traing_and_CV_error"),
#                                         manual_delt){
#   library(caret)
#   ##Pamr predict
#   #get train set and test set 
#   train_set <- train_and_test_dataset$train_set
#   test_set <- train_and_test_dataset$test_set
#   
#   if(is.null(rownames(train_set$x))){
#     rownames(train_set$x) <- paste0("gene_",1:nrow(train_set$x))
#   }
#   
#   if(is.null(rownames(test_set$x))){
#     rownames(test_set$x) <- paste0("gene_",1:nrow(test_set$x))
#   }
#   
#   #using training set to train model
#   training_model <- pamr.train(train_set,n.threshold = n_threshold)
#   
#   #cross_validation <- pamr.cv(fit = training_model,data = train_set,nfold = 10)
#   cross_validation <- pamr.cv_with_Cor(fit = training_model,data = train_set,nfold = 10)
#   # 
#   # pamr.plotcv(cross_validation)
#   
#   #determine delt value of pam
#   #you can choose manully or determine it by choosing the mean of the lowest prediction error rate and C-V error rate
#   if(auto_pam_delt){
#     #obtain the delt which give the lowest training error
#     errors <- setNames(training_model$errors,1:length(training_model$errors))
#     errors <- errors[length(training_model$errors):1]
#     min_error_from_last <- as.integer(names(which.min(errors)))
#     min_error_delt <- training_model$threshold[min_error_from_last]
#     
#     #obtain the delt which give the lowest training error
#     CV_errors <- setNames(cross_validation$error,1:length(cross_validation$error))
#     CV_errors <- CV_errors[length(cross_validation$error):1]
#     min_CV_error_from_last <- as.integer(names(which.min(CV_errors)))
#     min_CV_error_delt <- cross_validation$threshold[min_CV_error_from_last]
#     
#     switch(auto_delt_by,
#            min_training_error={
#              pam_delt <- min_error_delt
#            },
#            min_CV_error={
#              pam_delt <- min_CV_error_delt
#            },
#            mean_of_min_traing_and_CV_error={
#              pam_delt <- mean(c(min_error_delt,min_CV_error_delt))
#            },
#            stop("auto_delt_by arg must specified as one of (min_training_error min_CV_error mean_of_min_traing_and_CV_error)"))
#   }else{
#     if(missing(manual_delt)){
#       print(training_model)
#       print(cross_validation)
#       pam_delt <- as.numeric(readline(prompt = "which delt you choose? \n"))
#     }else{
#       message(paste0("You are manually specifying delt: ",manual_delt))
#       pam_delt <- manual_delt
#     }
#   }
#   test_predicted <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt)
#   train_predicted <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt)
#   
#   test_posterior <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt,type = "posterior")
#   train_posterior <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "posterior")
#   
#   ##correlation predicted
#   #assign geneid and genenames for training model,
#   #otherwise the pamr.listgenes would would throw error
#   train_set$geneid <- 1:nrow(train_set$x) #make the $geneid components in train_set the same as training_model$gene.subset
#   #if(!is.null(rownames(train_set$x))){
#   train_set$genenames <- rownames(train_set$x) #make the $genenames components in train_set the same as training_model$gene.subset
#   #}else{
#   #  train_set$genename <- as.character(1:nrow(train_set$x))
#   #}
#   
#   #pamr centroids
#   pamr_centroids <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "centroid")
#   pam_info <- capture.output(pamr_survived_genes <- pamr.listgenes(fit = training_model,
#                                                                    data = train_set,threshold = pam_delt,
#                                                                    genenames = TRUE)[,"name"])
#   
#   #if you want to show the info from pamr.listgenes, uncomment the command below
#   #cat(paste0(pam_info,collapse = "\n"))
#   
#   
#   # pamr_survived_genes <- pamr.listgenes(fit = training_model,data = train_set,threshold = pam_delt,genenames = TRUE)[,"id"]
#   if(length(pamr_survived_genes)<=1){
#     message("Number of survived genes less or equal to 1 produced which would lead to error in the following analysis
#             because structure of resulting data were changed!")
#     return(list(delt=pam_delt,gene_numbers=length(pamr_survived_genes),traning_model=training_model,CR_model=cross_validation,
#                 expr=train_and_test_dataset))
#   }
#   pamr_centroids <- pamr_centroids[pamr_survived_genes,]
#   
#   #correlation
#   pamr_centroids_test_cor <- cor(x = test_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
#   pamr_centroids_train_cor <- cor(x = train_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
#   
#   test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
#   train_cor_predicted <- apply(pamr_centroids_train_cor,1,function(x)names(which.max(x)))
#   
#   predict_res<-list(train_set=data.frame(pam_predict=train_predicted,
#                                          cor_predict=factor(train_cor_predicted,levels = levels(train_predicted)),
#                                          stringsAsFactors = FALSE),#set both column of data.frame as factor, so that
#                     #there would be the same columns in the table() function output, 
#                     #otherwise there would be an error in caret::confusionMatrix
#                     #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
#                     #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
#                     #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
#                     test_set=data.frame(pam_predict=test_predicted,
#                                         cor_predict=factor(test_cor_predicted,levels = levels(test_predicted)),
#                                         stringsAsFactors = FALSE)#set both column of data.frame as factor, so that
#                     #there would be the same columns in the table() function output, 
#                     #otherwise there would be an error in caret::confusionMatrix
#                     #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
#                     #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
#                     #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
#   )
#   kappas <- sapply(predict_res,function(x){
#     confusionMatrix(table(x))$overall["Kappa"]
#   })
#   
#   predict_res_error <- list(train_set=c(sum(train_predicted==train_and_test_dataset$train_set$y)/length(train_predicted),
#                                         sum(train_cor_predicted==train_and_test_dataset$train_set$y)/length(train_cor_predicted)),
#                             test_set=c(sum(test_predicted==train_and_test_dataset$test_set$y)/length(test_predicted),
#                                        sum(test_cor_predicted==train_and_test_dataset$test_set$y)/length(test_cor_predicted)))
#   predict_res_error$train_set <- setNames(predict_res_error$train_set,c("PAM_predicted","Cor_predicted"))
#   predict_res_error$test_set <- setNames(predict_res_error$test_set,c("PAM_predicted","Cor_predicted"))
#   
#   if(!(all(train_and_test_dataset$train_set$y%in%levels(train_predicted))|
#        all(train_and_test_dataset$test_set$y%in%levels(test_predicted))))
#     stop("The levels of train or test predicted results is not match the real class")
#   predict_res$train_set$real_class <- factor(train_and_test_dataset$train_set$y,levels = levels(train_predicted))
#   predict_res$test_set$real_class <- factor(train_and_test_dataset$test_set$y,levels = levels(test_predicted))
#   
#   confusion_matrices <- list(train_set=list(PAM_and_real=confusionMatrix(table(predict_res$train_set[,c("pam_predict",
#                                                                                                         "real_class")])),
#                                             Cor_and_real=confusionMatrix(table(predict_res$train_set[,c("cor_predict",
#                                                                                                         "real_class")]))),
#                              test_set=list(PAM_and_real=confusionMatrix(table(predict_res$test_set[,c("pam_predict",
#                                                                                                       "real_class")])),
#                                            Cor_and_real=confusionMatrix(table(predict_res$test_set[,c("cor_predict",
#                                                                                                       "real_class")]))))
#   
#   
#   # list(kappas=kappas,predict_res=predict_res,train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor,
#   #      centroid=pamr_centroids,delt=pam_delt,traning_and_testing=train_and_test_dataset)
#   
#   list(kappas=kappas,predict_res=predict_res,
#        pamr_posterior=list(train_set=train_posterior,test_set=test_posterior),
#        sample_cor=list(train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor),
#        centroid=pamr_centroids,delt=pam_delt,traning_and_testing=train_and_test_dataset,
#        error_rate=predict_res_error,predicted_confusions=confusion_matrices,
#        CV_res=cross_validation)
# }
# 

compare_pam_and_correlation <- function(train_and_test_dataset,auto_pam_delt=FALSE,n_threshold=100,
                                        auto_delt_by=c("min_training_error","min_CV_error","mean_of_min_traing_and_CV_error"),
                                        manual_delt,pamr_prior_class=c("euqal","class")){
  library(caret)
  ##Pamr predict
  #get train set and test set 
  train_set <- train_and_test_dataset$train_set
  test_set <- train_and_test_dataset$test_set
  
  if(is.null(rownames(train_set$x))){
    rownames(train_set$x) <- paste0("gene_",1:nrow(train_set$x))
  }
  
  if(is.null(rownames(test_set$x))){
    rownames(test_set$x) <- paste0("gene_",1:nrow(test_set$x))
  }
  
  #determine pamr_prior
  pamr_prior <- switch(pamr_prior_class,equal={
    prior <- 1/length(unique(train_set$y))
    rep(prior,length(unique(train_set$y)))
  },class=NULL)
  
  print(pamr_prior)
    
  #using training set to train model
  training_model <- pamr.train(train_set,n.threshold = n_threshold,prior = pamr_prior)
  
  cross_validation <- pamr.cv(fit = training_model,data = train_set,nfold = 10)
  #cross_validation <- pamr.cv_with_Cor(fit = training_model,data = train_set,nfold = 10) #the prior of pamr.cv is taken from 
                                                                                         #the results of pamr.train, 
                                                                                         #so we shouldn't assign prior here again
  # 
  # pamr.plotcv(cross_validation)
  
  #determine delt value of pam
  #you can choose manully or determine it by choosing the mean of the lowest prediction error rate and C-V error rate
  if(auto_pam_delt){
    #obtain the delt which give the lowest training error
    errors <- setNames(training_model$errors,1:length(training_model$errors))
    errors <- errors[length(training_model$errors):1]
    min_error_from_last <- as.integer(names(which.min(errors)))
    min_error_delt <- training_model$threshold[min_error_from_last]
    
    #obtain the delt which give the lowest training error
    CV_errors <- setNames(cross_validation$error,1:length(cross_validation$error))
    CV_errors <- CV_errors[length(cross_validation$error):1]
    min_CV_error_from_last <- as.integer(names(which.min(CV_errors)))
    min_CV_error_delt <- cross_validation$threshold[min_CV_error_from_last]
    
    switch(auto_delt_by,
           min_training_error={
             pam_delt <- min_error_delt
           },
           min_CV_error={
             pam_delt <- min_CV_error_delt
           },
           mean_of_min_traing_and_CV_error={
             pam_delt <- mean(c(min_error_delt,min_CV_error_delt))
           },
           stop("auto_delt_by arg must specified as one of (min_training_error min_CV_error mean_of_min_traing_and_CV_error)"))
  }else{
    if(missing(manual_delt)){
      print(training_model)
      print(cross_validation)
      pam_delt <- as.numeric(readline(prompt = "which delt you choose? \n"))
    }else{
      message(paste0("You are manually specifying delt: ",manual_delt))
      pam_delt <- manual_delt
    }
  }
  test_predicted <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt)
  train_predicted <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt)
  
  test_posterior <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt,type = "posterior")
  train_posterior <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "posterior")
  
  ##correlation predicted
  #assign geneid and genenames for training model,
  #otherwise the pamr.listgenes would would throw error
  train_set$geneid <- 1:nrow(train_set$x) #make the $geneid components in train_set the same as training_model$gene.subset
  #if(!is.null(rownames(train_set$x))){
  train_set$genenames <- rownames(train_set$x) #make the $genenames components in train_set the same as training_model$gene.subset
  #}else{
  #  train_set$genename <- as.character(1:nrow(train_set$x))
  #}
  
  #pamr centroids
  pamr_centroids <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "centroid")
  pam_info <- capture.output(pamr_survived_genes <- pamr.listgenes(fit = training_model,
                                                                   data = train_set,threshold = pam_delt,
                                                                   genenames = TRUE)[,"name"])
  
  #if you want to show the info from pamr.listgenes, uncomment the command below
  #cat(paste0(pam_info,collapse = "\n"))
  
  
  # pamr_survived_genes <- pamr.listgenes(fit = training_model,data = train_set,threshold = pam_delt,genenames = TRUE)[,"id"]
  if(length(pamr_survived_genes)<=1){
    message("Number of survived genes less or equal to 1 produced which would lead to error in the following analysis
            because structure of resulting data were changed!")
    return(list(delt=pam_delt,gene_numbers=length(pamr_survived_genes),traning_model=training_model,CR_model=cross_validation,
                expr=train_and_test_dataset))
  }
  pamr_centroids <- pamr_centroids[pamr_survived_genes,]
  
  #correlation
  pamr_centroids_test_cor <- cor(x = test_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
  pamr_centroids_train_cor <- cor(x = train_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
  
  test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
  train_cor_predicted <- apply(pamr_centroids_train_cor,1,function(x)names(which.max(x)))
  
  predict_res<-list(train_set=data.frame(pam_predict=train_predicted,
                                         cor_predict=factor(train_cor_predicted,levels = levels(train_predicted)),
                                         stringsAsFactors = FALSE),#set both column of data.frame as factor, so that
                    #there would be the same columns in the table() function output, 
                    #otherwise there would be an error in caret::confusionMatrix
                    #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
                    #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
                    #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
                    test_set=data.frame(pam_predict=test_predicted,
                                        cor_predict=factor(test_cor_predicted,levels = levels(test_predicted)),
                                        stringsAsFactors = FALSE)#set both column of data.frame as factor, so that
                    #there would be the same columns in the table() function output, 
                    #otherwise there would be an error in caret::confusionMatrix
                    #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
                    #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
                    #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
  )
  kappas <- sapply(predict_res,function(x){
    confusionMatrix(table(x))$overall["Kappa"]
  })
  
  predict_res_error <- list(train_set=c(sum(train_predicted!=train_and_test_dataset$train_set$y)/length(train_predicted),
                                        sum(train_cor_predicted!=train_and_test_dataset$train_set$y)/length(train_cor_predicted)),
                            test_set=c(sum(test_predicted!=train_and_test_dataset$test_set$y)/length(test_predicted),
                                       sum(test_cor_predicted!=train_and_test_dataset$test_set$y)/length(test_cor_predicted)))
  predict_res_error$train_set <- setNames(predict_res_error$train_set,c("PAM_predicted","Cor_predicted"))
  predict_res_error$test_set <- setNames(predict_res_error$test_set,c("PAM_predicted","Cor_predicted"))
  
  if(!(all(train_and_test_dataset$train_set$y%in%levels(train_predicted))|
       all(train_and_test_dataset$test_set$y%in%levels(test_predicted))))
    stop("The levels of train or test predicted results is not match the real class")
  predict_res$train_set$real_class <- factor(train_and_test_dataset$train_set$y,levels = levels(train_predicted))
  predict_res$test_set$real_class <- factor(train_and_test_dataset$test_set$y,levels = levels(test_predicted))
  
  confusion_matrices <- list(train_set=list(PAM_and_real=confusionMatrix(table(predict_res$train_set[,c("pam_predict",
                                                                                                        "real_class")])),
                                            Cor_and_real=confusionMatrix(table(predict_res$train_set[,c("cor_predict",
                                                                                                        "real_class")]))),
                             test_set=list(PAM_and_real=confusionMatrix(table(predict_res$test_set[,c("pam_predict",
                                                                                                      "real_class")])),
                                           Cor_and_real=confusionMatrix(table(predict_res$test_set[,c("cor_predict",
                                                                                                      "real_class")]))))
  
  
  # list(kappas=kappas,predict_res=predict_res,train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor,
  #      centroid=pamr_centroids,delt=pam_delt,traning_and_testing=train_and_test_dataset)
  
  list(kappas=kappas,predict_res=predict_res,
       pamr_posterior=list(train_set=train_posterior,test_set=test_posterior),
       sample_cor=list(train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor),
       centroid=pamr_centroids,delt=pam_delt,traning_and_testing=train_and_test_dataset,
       error_rate=predict_res_error,predicted_confusions=confusion_matrices,
       CV_res=cross_validation,training_model=training_model)
}




n_times_compare <- function(expr_with_label,n_times,train_proportion,auto_pam_delt,n_threshold,prior,...){
  compare_res <- list()
  i <- 1
  while(i <= n_times){
    message(paste0("The ",i,"th compare!"))
    train_and_test_dat <- produce_train_test_set(expr_with_label = expr_with_label,train_proportion = train_proportion)
    # test_error <- PAM_test_error(train_and_test_dataset = train_and_test_dat,
    #                              auto_pam_delt = auto_pam_delt,n_threshold = n_threshold)
    com_res <- compare_pam_and_correlation(train_and_test_dataset = train_and_test_dat,auto_pam_delt = auto_pam_delt,
                                           n_threshold = n_threshold,pamr_prior_class = prior,...)
    compare_res <- c(compare_res,list(com_res))
    i <- i+1
  }
  setNames(compare_res,paste0("The ",1:n_times,"st errors"))
}


gene_number_and_kappa <- function(compare_res,supposed_compare_res_length=7,title){
  #find failed samples
  component_length <- sapply(compare_res,length)
  compare_res <- compare_res[(component_length==supposed_compare_res_length)]
  all_kappas <- t(sapply(compare_res,function(x)x$kappas))
  all_cen_numbers <- log10(sapply(compare_res,function(x)nrow(x$centroid)))
  # plot(all_cen_numbers,all_kappas[,1],type="l",col="red",lwd=2,xlab="Number of Centroid Gene",ylab="Kappa",font.lab=2,cex=1.5,
  #      ylim=c(min(all_kappas),max(all_kappas)))
  # lines(x = all_cen_numbers,all_kappas[,2],type = "l",col="blue",lwd=2)
  # legend("bottomright",legend = c("Traning","Testing"),col = c("red","blue"),lwd = 2,cex = 1.5,bty="n")
  
  plot(all_cen_numbers,all_kappas[,1],col=ggplot2::alpha("red",alpha = 0.3),xlab="Number of Centroid Gene (log10)",
       ylab="Kappa",font.lab=2,cex=1.5,pch=16,
       ylim=c(min(all_kappas),max(all_kappas)),main=title)
  points(x = all_cen_numbers,all_kappas[,2],col=ggplot2::alpha("blue",alpha = 0.3),cex=1.5,pch=16)
  legend("bottomright",legend = c("Traning","Testing"),
         col = c(ggplot2::alpha("red",alpha = 0.3),ggplot2::alpha("blue",alpha = 0.3)),cex = 1.5,bty="n",pch = 16)
}



#****************************************************************
#******Compare ClaNC classifying and pAM + spearman results******
#****************************************************************

ClaNC_clasification_res <- function(train_and_test_dataset,class_label,prior="class",
                                    CV_gene_number=1:50,auto_active_genes=TRUE){
  #get train set and test set 
  train_set <- train_and_test_dataset$train_set
  test_set <- train_and_test_dataset$test_set
  
  if(is.null(rownames(train_set$x))){
    rownames(train_set$x) <- paste0("gene_",1:nrow(train_set$x))
  }
  if(is.null(rownames(test_set$x))){
    rownames(test_set$x) <- paste0("gene_",1:nrow(test_set$x))
  }
  
  build_model <- ClaNC_trained_and_builed_model(expr_set = train_set,gene_names = rownames(train_set$x),
                                               sample_class = train_set$y,class_label = class_label,
                                               prior = prior,CV_gene_number = CV_gene_number,
                                               auto_active_genes = auto_active_genes)
  
  Training_predicted_res <- ClaNC_predicted(ClaNC_builded = build_model$builded,expr_set = train_set,
                                            gene_names = rownames(train_set$x),class_label = class_label)
  Testing_predicted_res <- ClaNC_predicted(ClaNC_builded = build_model$builded,expr_set = test_set,
                                           gene_names = rownames(test_set$x),class_label = class_label)
  train_errors <- testClanc(data = train_set$x,id = train_set$y,geneNames = rownames(train_set$x),fit = build_model$builded)
  test_errors <- testClanc(data = test_set$x,id = test_set$y,geneNames = rownames(test_set$x),fit = build_model$builded)
  
  
  list(training_predict_res=Training_predicted_res,testing_predicted_res=Testing_predicted_res,
       gene_fetures=rownames(build_model$builded$cntrds),predicting_model=build_model,
       errors=list(train_error_number=setNames(train_errors,class_label),
                   test_error_number=setNames(test_errors,class_label)),
	   class_label=class_label)
}

ClaNC_trained_and_builed_model <- function(expr_set,gene_names,sample_class,class_label,prior="class",
                                           CV_gene_number=1:50,auto_active_genes=TRUE){
  if(all(order(CV_gene_number)!=(1:length(CV_gene_number))))
    stop("CV_gene_number must be given as increasing order, such as 1, 2, 3, 4, 5, 6, ...")
  expr_dat <- as.matrix(expr_set$x)
  
  #Cross validation
  CV_res <- cvClanc(data = expr_dat,id = sample_class,prior = prior,active = CV_gene_number)#default up to 5 folds CV
  if(auto_active_genes){
    message(paste0("You specified choosing active gene automatically,\n",
                   "which means the number of active genes would be choosed as the smallest one to give minimum CV error"))
    
    #this command make sure we use as least as genes which give
    #the minimum CV errors
    active_gene_number <- CV_gene_number[which.min(CV_res$overallErrors)[1]]
  }else{
    plot(CV_gene_number,CV_res$overallErrors,type = "l", lwd = 2, col = "blue", xlab = 
           "Number of features", ylab = "CV error")
    active_gene_number <- readline(prompt = "Which number of active gene you choose: ")
    active_gene_number <- as.integer(active_gene_number)
  }
  
  #build the prediction model
  trained_model <- trainClanc(data = expr_dat,id = sample_class,geneNames = gene_names,prior = prior)
  builded_model <- buildClanc(data = expr_dat,id = sample_class,cNames = class_label,
                              train = trained_model,active = active_gene_number)
  
  list(trained=trained_model,builded=builded_model,active_gene_number=active_gene_number,CV_res=CV_res)
}

ClaNC_predicted <- function(ClaNC_builded,expr_set,gene_names,class_label){
  predict_res <- predictClanc(data = expr_set$x,geneNames = gene_names,fit = ClaNC_builded)
  
  #predicted_label <- as.character(factor(predict_res,levels = unique(predict_res), labels = class_label))
  predicted_label <- class_label[predict_res]
  
  setNames(predict_res,predicted_label)
}


# #a <- produce_train_test_set(expr_with_label = Pomeroy_datset,train_proportion = 2/3)
# setwd("d:/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
# source("clanc.R")
# Y_train = read.delim("data/trainExample.txt", header = T)
# Y_test = read.delim("data/testExample.txt", header = T)
# rownames(Y_train) = Y_train[, 1]
# rownames(Y_test) = Y_test[, 1]
# Y_train = as.matrix(Y_train[, -1])
# Y_test = as.matrix(Y_test[, -1])
# 
# dat_predict <- list(train_set=list(x=Y_train,y=rep(1:4, each = 10)),
#                     test_set=list(x=Y_test,y=rep(1:4, each = 10)))
# 
# predicted_res <- ClaNC_clasification_res(train_and_test_dataset = dat_predict,class_label = letters[1:4],prior = "class",
#                                          CV_gene_number = 1:50,auto_active_genes = TRUE)




# compare_ClaNC_and_PAM50 <- function(expr_with_label, train_proportion = 2/3, 
#                                     prior="class", CV_gene_number = 1:50, auto_active_genes = TRUE,...){
compare_ClaNC_and_PAM50 <- function(expr_with_label, train_proportion = 2/3, already_train_and_test_set=FALSE,
                                    prior="class", CV_gene_number = 1:50, auto_active_genes = TRUE,...){
  
  if(already_train_and_test_set){
    train_and_test <- expr_with_label
    class_label <- levels(as.factor(train_and_test$train_set$y))
    
    #change $y in train_and_test$train_set and $test_set to numeric to satisfy ClaNC functions
    train_and_test <- lapply(train_and_test, function(x) {
      y <- as.integer(factor(x$y,levels = sort(unique(x$y))))
      list(x=x$x,y=y)
    })
    
  }else{
    #determine the class labels
    expr_with_numeric_label <- list(x=expr_with_label$x,y=as.integer(as.factor(expr_with_label$y)))
    
    
    train_and_test <- produce_train_test_set(expr_with_label = expr_with_numeric_label,train_proportion = train_proportion)
    
    class_label <- levels(as.factor(expr_with_label$y)) #The class_label is corresponded with the label order of $y in train_and_test
    #and it is not corresponded with expr_with_label
    #This is right since in the function produce_train_test_set, the 
    #the label order is changed and it is consistent with expr_with_label
    
    #change $y in train_and_test$train_set and $test_set to numeric to satisfy ClaNC functions
    train_and_test <- lapply(train_and_test, function(x) {
      y <- as.integer(factor(x$y,levels = unique(x$y))) 
      ##the command above and below would obtain the same result in this situation
      #y <- as.integer(factor(x$y,levels = sort(unique(x$y))))
      list(x=x$x,y=y)
    })
  }
  
  
  
  # message("The order of class label must match the order of unique(expr_with_label) !!!")
  
  ClaNC_predict_res <- ClaNC_clasification_res(train_and_test_dataset = train_and_test,class_label = class_label,
                                               prior = prior,CV_gene_number = CV_gene_number,
                                               auto_active_genes = auto_active_genes)
  
  #use the ClaNC selected genes
  train_and_test_for_PAM <- lapply(train_and_test,function(x,genes,class_label){
    selected_gene_dat <- x$x[genes,]
    y <- class_label[x$y]
    list(x=selected_gene_dat,y=y)
  },genes=ClaNC_predict_res$gene_fetures,class_label=class_label)
  
  #PAM_and_spearman_res <- compare_pam_and_correlation(train_and_test_dataset = train_and_test_for_PAM, ...)
  
  PAM_and_spearman_res <- compare_pam_and_correlation_for_cor(train_and_test_dataset = train_and_test_for_PAM,
                                                              pamr_prior_class = prior,...)
  
  ##train_set kappa
  #ClaNC VS centroids + cor
  ClaNC_and_PAM_plus_spearman_train <- data.frame(ClaNC=factor(names(ClaNC_predict_res$training_predict_res),
                                                               levels = unique(class_label)),
                                                  PAM_and_Spearman=factor(PAM_and_spearman_res$predict_res$train_set$cor_predict,
                                                                          levels = unique(class_label)))
  #ClaNC + PAM predict function VS centroids + cor
  # ClaNC_PAM_and_PAM_plus_spearman_train <- PAM_and_spearman_res$train_set[,c("pam_predict","cor_predict")]
  # colnames(ClaNC_PAM_and_PAM_plus_spearman_train) <- c("ClaNC","PAM_and_Spearman")
  
  ##test_set kappa
  #ClaNC VS centroids + cor
  ClaNC_and_PAM_plus_spearman_test <- data.frame(ClaNC=factor(names(ClaNC_predict_res$testing_predicted_res),
                                                              levels = unique(class_label)),
                                                 PAM_and_Spearman=factor(PAM_and_spearman_res$predict_res$test_set$cor_predict,
                                                                         levels = unique(class_label)))
  #ClaNC + PAM predict function VS centroids + cor
  # ClaNC_PAM_and_PAM_plus_spearman_test <- PAM_and_spearman_res$test_set[,c("pam_predict","cor_predict")]
  # colnames(ClaNC_PAM_and_PAM_plus_spearman_test) <- c("ClaNC","PAM_and_Spearman")
  
  
  train_kappa <- confusionMatrix(table(ClaNC_and_PAM_plus_spearman_train))$overall["Kappa"]
  test_kappa <- confusionMatrix(table(ClaNC_and_PAM_plus_spearman_test))$overall["Kappa"]
  
  consistency_between_ClaNC_and_mean_plus_cor <- c(train_kappa,test_kappa)
  
  ClaNC_predicted_and_real_class <- list(train_set=data.frame(ClaNC=factor(names(ClaNC_predict_res$training_predict_res),
                                                                           levels = unique(class_label)),
                                                              real_class=factor(train_and_test_for_PAM$train_set$y,
                                                                                levels = unique(class_label))),
                                         test_set=data.frame(ClaNC=factor(names(ClaNC_predict_res$testing_predicted_res),
                                                                          levels = unique(class_label)),
                                                             real_class=factor(train_and_test_for_PAM$test_set$y,
                                                                               levels = unique(class_label))))
  predicted_confusion <- lapply(ClaNC_predicted_and_real_class,function(x){
    confusionMatrix(table(x))
  })
  
  
  # predicted_confusions <- list(train_set=list(ClaNC_and_PAM_plus_cor=confusionMatrix(table(ClaNC_and_PAM_plus_spearman_train)),
  #                                             ClaNC_PAM_and_PAM_plus_cor=),
  #                              test_set=list(ClaNC_and_PAM_plus_cor=confusionMatrix(ClaNC_and_PAM_plus_spearman_test),
  #                                            ClaNC_PAM_and_PAM_plus_cor=))
  
  ClaNC_predict_res$predicted_confusions <- predicted_confusion
  
  list(ClaNC_and_mean_plus_cor=setNames(consistency_between_ClaNC_and_mean_plus_cor,c("train_set.Kappa","test_set.Kappa")),
       ClaNC_plus_mean_PAM_and_mean_plus_cor=PAM_and_spearman_res$kappas,
       ClaNC_res=ClaNC_predict_res,PAM_res=PAM_and_spearman_res)
}

n_times_CLaNC_and_PAM_com <- function(expr_with_label,ntimes=50,show_message=TRUE,...){
  
  compare_res <- list()
  i <- 1
  message_function <- ifelse(show_message,base::get("("),suppressMessages)
  
  cat_suffix <- ifelse(show_message,"\n","")
  cat(paste0("\r",format(paste0("Total ",ntimes," times: "),width = 17,justify = "r"),
             format(paste0(rep("=",floor(((i-1)/ntimes)*20),collapse=""),collapse = ""),width = 20,justify = "l"),
             ": ",sprintf("%5.1f",(i-1)*100/ntimes),"% finished!",cat_suffix,sep = ""))
  
  while(i <= ntimes){
    #message(paste0("The ",i,"th compare!"))
    message_function(print_info <- capture.output(ith_compare <- compare_ClaNC_and_PAM50(expr_with_label = expr_with_label,...)))
    ##To see the printing infomation, uncomment the command above and use the command below
    #ith_compare <- compare_ClaNC_and_PAM50(expr_with_label = expr_with_label,...)
    
    compare_res <- c(compare_res,list(ith_compare))
    i <- i + 1
    
    # message(paste0(format(paste0("Total ",ntimes," times: "),width = 17,justify = "r"),
    #                format(paste0(rep("=",floor(((i-1)/ntimes)*20),collapse=""),collapse = ""),width = 20,justify = "l"),
    #                ": ",sprintf("%5.1f",(i-1)*100/ntimes),"% finished!\n",sep = ""))
    cat(paste0("\r",format(paste0("Total ",ntimes," times: "),width = 17,justify = "r"),
               format(paste0(rep("=",floor(((i-1)/ntimes)*20),collapse=""),collapse = ""),width = 20,justify = "l"),
               ": ",sprintf("%5.1f",(i-1)*100/ntimes),"% finished!",cat_suffix,sep = ""))
  }
  cat("\n")
  setNames(compare_res,paste0("The ",1:ntimes,"st errors"))
}
# setwd("d:/Shi_lab/Breast_cancer/ClaNC_classifier/clanc_share/")
# source("clanc.R")

cvClanc <- function(data, id, prior = "equal", active = 1:10, gui = F, prntFrm = NULL) {
  ## data: expression data.  (m x n) matrix of class numeric.
  ## id: class id's.  n vector of class numeric.
  ## prior: either a vector of length p, "equal", or "class".  if "equal", then equal probabilities 
  ##   will be used for each class.  if "class", then the proportions of each class in the 
  ##   training data will be used as the prior.
  ## active: how many active features to consider?  can either be a single number or a vector 
  ##   containing a range of values to consider.
  ## gui: indicates whether call is coming from gui.  typical user should leave default.
  ## prntFrm: frame in which to print if called from gui.  typical user should leave default.
  
  cvIdx = balancedFolds(id, 5)
  
  m = nrow(data)
  n = ncol(data)
  p = length(unique(id))
  nn = as.numeric(table(id))
  folds = length(cvIdx)
  
  if(is.numeric(prior)) {
    if(length(prior) != p | sum(prior) != 1)
      stop("Invalid prior.")
    pi.k = prior
  } else {
    if(prior == "equal")
      pi.k = rep(1 / p, p)
    else if(prior == "class")
      pi.k = nn / n
  }
  
  if(is.matrix(active)) {
    d = nrow(active)
  } else {
    d = length(active)
  }
  
  if(gui) {
    printClanc = function(msg) {
      assign("shareMsg", msg, prntFrm)
      eval(expression(postMsg(shareMsg)), prntFrm)
    }
  } else {
    printClanc = cat
  }
  
  cv.error = array(rep(0, d * folds * p), dim = c(d, folds, p))
  
  cv.err.cnt.cls = matrix(NA, nrow = d, ncol = p)
  cv.err.prpn.cls = matrix(NA, nrow = d, ncol = p)
  n.features.cls = matrix(NA, nrow = d, ncol = p)
  n.features.ttl = rep(NA, d)
  predicted_class <- matrix(0,nrow = length(id),ncol = d)
  colnames(predicted_class) <- paste0(1:d,"_active_genes")
  rownames(predicted_class) <- colnames(data)
  
  ID = model.matrix(~ factor(id) - 1)
  dimnames(ID) = list(NULL, names(nn))
  
  ## cross validation
  printClanc("CV:")
  for(i in 1:folds) {
    printClanc(i)
    
    ## form initial statistics
    v = length(cvIdx[[i]])
    X = data[, -cvIdx[[i]]]
    Y = data[, cvIdx[[i]]]
    jd = id[-cvIdx[[i]]]
    truth = id[cvIdx[[i]]]
    JD = ID[-cvIdx[[i]], ]
    
    mm = table(jd)
    m.k = sqrt(1 / mm - 1 / (n - v))
    
    ## pooled standard deviations
    p.sd = pooledSDClanc(X, JD)
    
    ## class- and overall-centroids
    cntrd.k = scale(X %*% JD, center = F, scale = mm)
    cntrd.o = drop(X %*% rep(1 / (n - v), n - v))
    
    ## form statistics and order them
    d.k = scale((cntrd.k - cntrd.o) / p.sd, center = F, scale = m.k)
    d.k.ord = orderStatsClanc(d.k = d.k)
    
    ## select genes, update inactive centroid components
    for(j in 1:d) {
      if(is.matrix(active))
        aa = active[j, ]
      else
        aa = active[j]
      
      selected = selectClanc(d.k = d.k, d.k.ord = d.k.ord, active = aa)
      active.idx = (1:m)[drop(selected %*% rep(1, p)) != 0]
      
      cntrds = cntrd.k[active.idx, ]
      for(k in 1:p)
        cntrds[selected[active.idx, k] == 0, k] = cntrd.o[active.idx][selected[active.idx, k] == 0]
      
      ## classify test sample and assess error status
      for(k in 1:v) {
        dd = distClanc(data = Y[active.idx, k], cntrds = cntrds, sd = p.sd[active.idx], prior = pi.k)
        
        if(match(min(dd), dd) != truth[k])
          cv.error[j, i, truth[k]] = cv.error[j, i, truth[k]] + 1
        
        predicted_class[cvIdx[[i]][k],j] <- which.min(dd)
      }
    }
  }
  
  ## record numbers and proportions of errors
  for(i in 1:p) {
    if(d > 1)
      cv.err.cnt.cls[, i] = apply(cv.error[, , i], 1, sum)
    else
      cv.err.cnt.cls[, i] = sum(cv.error[, , i])
    cv.err.prpn.cls[, i] = cv.err.cnt.cls[, i] / nn[i]
  }
  cv.err.cnt.ttl = apply(cv.err.cnt.cls, 1, sum)
  cv.err.prpn.ttl = cv.err.cnt.ttl / n
  
  printClanc("\n")
  return(list("classErrors" = cv.err.prpn.cls, "overallErrors" = cv.err.prpn.ttl, "prior" = pi.k,
              CV_predict_matrix=predicted_class,sample_real_class=id))
}

"balancedFolds" <- 
function(y, nfolds = 5) {
  permuteRows = function(x) {
    dd = dim(x)
    n = dd[1]
    p = dd[2]
    mm = runif(length(x)) + rep(seq(n) * 10, rep(p, n))
    matrix(t(x)[order(mm)], n, p, byrow = T)
  }

  totals = table(y)
  fmax = max(totals)
  nfolds = min(nfolds, fmax)
  folds = as.list(seq(nfolds))
  yids = split(seq(y), y)
  bigmat = matrix(NA, ceiling(fmax / nfolds) * nfolds, length(totals))
  for(i in seq(totals))
    bigmat[seq(totals[i]), i] = sample(yids[[i]])
  smallmat = matrix(bigmat, nrow = nfolds)
  smallmat = permuteRows(t(smallmat))
  res = vector("list", nfolds)
  for(j in 1:nfolds) {
    jj = !is.na(smallmat[, j])
    res[[j]] = smallmat[jj, j]
  }

  return(res)
}


ClaNC_validation_MCC <- function(ClaNC_res){
  active_genes <- ClaNC_res$predicting_model$active_gene_number
  predicted_class <- ClaNC_res$predicting_model$CV_res$CV_predict_matrix[,paste0(active_genes,"_active_genes")]
  real_class <- ClaNC_res$predicting_model$CV_res$sample_real_class
  mltools::mcc(preds = predicted_class,actuals = real_class)
}

PAM_validation_MCC <- function(PAM_res){
  delt <- PAM_res$delt
  predicted_class <- PAM_res$CV_res$yhat[,PAM_res$CV_res$threshold==delt]
  real_class <- PAM_res$predict_res$train_set$real_class
  mltools::mcc(preds = predicted_class,actuals = real_class)
}

Cor_validation_MCC <- function(PAM_res){
  delt <- PAM_res$delt
  predicted_class <- PAM_res$CV_res$Cor_yhat[,PAM_res$CV_res$threshold==delt]
  real_class <- PAM_res$predict_res$train_set$real_class
  mltools::mcc(preds = predicted_class,actuals = real_class)
}


pamr.cv_with_Cor <- function (fit, data, nfold = NULL, folds = NULL, ...) 
{
  x <- data$x[fit$gene.subset, fit$sample.subset]
  if (!is.null(data$y) & !is.null(data$proby)) {
    stop("Must have exactly one of y and  proby  present in the data object")
  }
  y <- NULL
  proby <- NULL
  if (!is.null(fit$y)) {
    y <- factor(fit$y[fit$sample.subset])
  }
  if (!is.null(fit$proby)) {
    proby <- fit$proby[fit$sample.subset, ]
  }
  this.call <- match.call()
  junk <- nsccv_with_cor(x, y = y, proby = proby, object = fit, nfold = nfold, 
                         folds = folds, survival.time = data$survival.time, censoring.status = data$censoring.status, 
                         ngroup.survival = fit$ngroup.survival, problem.type = fit$problem.type, 
                         ...)
  junk$call <- this.call
  junk$sample.subset <- fit$sample.subset
  class(junk) = "pamrcved"
  junk
}


nsccv_with_cor <- function (x, y = NULL, proby = NULL, nfold = min(table(y)), 
                            folds = NULL, threshold = NULL, threshold.scale = NULL, 
                            survival.time = NULL, censoring.status = NULL, ngroup.survival = NULL, 
                            prior, object, ...) 
{
  this.call <- match.call()
  argy <- y
  if (is.null(y)) {
    y <- as.factor(apply(proby, 1, which.is.max))
  }
  n <- length(y)
  if (is.null(nfold) & is.null(survival.time)) {
    nfold <- min(table(y))
  }
  if (is.null(nfold) & !is.null(survival.time)) {
    nfold <- 10
  }
  if (is.null(survival.time)) {
    if (is.null(folds)) {
      folds <- pamr:::balanced.folds(y)
    }
  }
  if (!is.null(survival.time)) {
    if (is.null(folds)) {
      folds <- split(sample(1:n), rep(1:nfold, length = n))
    }
  }
  nfold <- length(folds)
  if (missing(prior)) {
    if (missing(object)) 
      prior <- table(y)/n
    else prior <- object$prior
  }
  if (missing(threshold)) {
    if (missing(object)) 
      stop("Must either supply threshold argument, or an nsc object")
    else {
      threshold <- object$threshold
      threshold.scale <- object$threshold.scale
      se.scale <- object$se.scale
    }
  }
  n.threshold <- length(threshold)
  yhat <- rep(list(y), n.threshold)
  names(yhat) <- paste(seq(n.threshold))
  yhat <- data.frame(yhat)
  
  # Cor_yhat <- matrix(0,nrow = ncol(x),ncol = n.threshold)
  # rownames(Cor_yhat) <- colnames(x)
  # colnames(Cor_yhat) <-  1:n.threshold
  Cor_yhat <- rep(list(y), n.threshold)
  names(Cor_yhat) <-  paste(seq(n.threshold))
  Cor_yhat <-  data.frame(Cor_yhat)
  
  
  n.class <- table(y)
  prob <- array(1, c(n, length(n.class), n.threshold))
  size <- double(n.threshold)
  hetero <- object$hetero
  cv.objects = vector("list", nfold)
  for (ii in 1:nfold) {
    cat("Fold", ii, ":")
    a <- pamr:::nsc(x[, -folds[[ii]]], y = argy[-folds[[ii]]], 
                    x[, folds[[ii]], drop = FALSE], proby = proby[-folds[[ii]], 
                                                                  ], threshold = threshold, threshold.scale = threshold.scale, 
                    se.scale = se.scale, prior = prior, hetero = hetero, 
                    ..., remove.zeros = FALSE)
    size <- size + a$nonzero
    prob[folds[[ii]], , ] <- a$prob
    yhat[folds[[ii]], ] <- a$yhat
    cat("\n")
    cv.objects[[ii]] = a
    
    training_model <- pamr.train(data = list(x=x[, -folds[[ii]]],y=argy[-folds[[ii]]]),
                                 threshold = threshold)
    for(Nth_threshold in 1:length(threshold)){
      delt <- threshold[Nth_threshold]
      pamr_centroids <- pamr.predict(fit = training_model,newx = x[, folds[[ii]], drop = FALSE],threshold = delt,type = "centroid")
      
      pam_info <- capture.output(pamr_survived_genes <- pamr.listgenes(fit = training_model,
                                                                       data = list(x=x[, -folds[[ii]],drop=FALSE],y=argy[-folds[[ii]]],
                                                                                   geneid = 1:nrow(x),
                                                                                   genenames = rownames(x)),
                                                                       threshold = delt,
                                                                       genenames = TRUE)[,"name"])
      if(length(pamr_survived_genes)<=1){
        warning("Number of survived genes less or equal to 1, which lead to samples could not be predicted by Spearman Correlation!")
        Cor_yhat[folds[[ii]],Nth_threshold] <- NA
      }else{
        pamr_centroids <- pamr_centroids[pamr_survived_genes,]
        #correlation
        pamr_centroids_test_cor <- cor(x = x[pamr_survived_genes, folds[[ii]], drop = FALSE],y = pamr_centroids,method = "spearman")
        
        test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
        
        Cor_yhat[folds[[ii]],Nth_threshold] <- test_cor_predicted
      }
    }
    
  }
  if (missing(object)) 
    size <- round(size/nfold)
  else size <- object$nonzero
  error <- rep(NA, n.threshold)
  loglik <- error
  pvalue.survival <- error
  pvalue.survival.func <- function(group, survival.time, censoring.status, 
                                   ngroup.survival) {
    temp <- coxph(Surv(survival.time, censoring.status) ~ 
                    as.factor(group))
    loglik <- 2 * (temp$loglik[2] - temp$loglik[1])
    return(1 - pchisq(loglik, ngroup.survival - 1))
  }
  if (!is.null(proby)) {
    proby.temp <- proby
  }
  else if (!is.null(survival.time)) {
    proby.temp <- pamr.surv.to.class2(survival.time, censoring.status, 
                                      n.class = ngroup.survival)$prob
  }
  for (i in 1:n.threshold) {
    if (is.null(survival.time) & is.null(proby)) {
      error[i] <- sum(yhat[, i] != y)/n
    }
    if (!is.null(survival.time)) {
      temp <- c(yhat[, i], names(table(y)))
      Yhat <- model.matrix(~factor(temp) - 1, data = list(y = temp))
      Yhat <- Yhat[1:length(yhat[[ii]]), ]
      error[i] <- (length(yhat[, i]) - sum(Yhat * proby.temp))/n
    }
    if (is.null(survival.time)) {
      loglik[i] <- sum(log(prob[, , i][cbind(seq(1, n), 
                                             unclass(y))]))/n
    }
    if (!is.null(survival.time)) {
      pvalue.survival[i] <- pvalue.survival.func(yhat[, 
                                                      i], survival.time, censoring.status, ngroup.survival)
    }
  }
  obj <- list(threshold = threshold, error = error, loglik = loglik, 
              size = size, yhat = yhat, y = y, prob = prob, folds = folds, 
              cv.objects = cv.objects, pvalue.survival = pvalue.survival, 
              call = this.call,Cor_yhat=Cor_yhat)
  class(obj) <- "nsccv"
  obj
}

compare_pam_and_correlation_for_cor <- function(train_and_test_dataset,auto_pam_delt=FALSE,n_threshold=100,
                                                auto_delt_by=c("min_training_error","min_CV_error","mean_of_min_traing_and_CV_error"),
                                                manual_delt,pamr_prior_class=c("euqal","class")){
  library(caret)
  ##Pamr predict
  #get train set and test set 
  train_set <- train_and_test_dataset$train_set
  test_set <- train_and_test_dataset$test_set
  
  if(is.null(rownames(train_set$x))){
    rownames(train_set$x) <- paste0("gene_",1:nrow(train_set$x))
  }
  
  if(is.null(rownames(test_set$x))){
    rownames(test_set$x) <- paste0("gene_",1:nrow(test_set$x))
  }
  
  #determine pamr_prior
  pamr_prior <- switch(pamr_prior_class,equal={
    prior <- 1/length(unique(train_set$y))
    rep(prior,length(unique(train_set$y)))
  },class=NULL)
  
  print(pamr_prior)
  
  #using training set to train model
  training_model <- pamr.train(train_set,n.threshold = n_threshold,prior = pamr_prior)
  
  #cross_validation <- pamr.cv(fit = training_model,data = train_set,nfold = 10)
  cross_validation <- pamr.cv_with_Cor(fit = training_model,data = train_set,nfold = 10) #the prior of pamr.cv is taken from 
  #the results of pamr.train, 
  #so we shouldn't assign prior here again
  # 
  # pamr.plotcv(cross_validation)
  
  #determine delt value of pam
  #you can choose manully or determine it by choosing the mean of the lowest prediction error rate and C-V error rate
  if(auto_pam_delt){
    #obtain the delt which give the lowest training error
    errors <- setNames(training_model$errors,1:length(training_model$errors))
    errors <- errors[length(training_model$errors):1]
    min_error_from_last <- as.integer(names(which.min(errors)))
    min_error_delt <- training_model$threshold[min_error_from_last]
    
    #obtain the delt which give the lowest training error
    CV_errors <- setNames(cross_validation$error,1:length(cross_validation$error))
    CV_errors <- CV_errors[length(cross_validation$error):1]
    min_CV_error_from_last <- as.integer(names(which.min(CV_errors)))
    min_CV_error_delt <- cross_validation$threshold[min_CV_error_from_last]
    
    switch(auto_delt_by,
           min_training_error={
             pam_delt <- min_error_delt
           },
           min_CV_error={
             pam_delt <- min_CV_error_delt
           },
           mean_of_min_traing_and_CV_error={
             pam_delt <- mean(c(min_error_delt,min_CV_error_delt))
           },
           stop("auto_delt_by arg must specified as one of (min_training_error min_CV_error mean_of_min_traing_and_CV_error)"))
  }else{
    if(missing(manual_delt)){
      print(training_model)
      print(cross_validation)
      pam_delt <- as.numeric(readline(prompt = "which delt you choose? \n"))
    }else{
      message(paste0("You are manually specifying delt: ",manual_delt))
      pam_delt <- manual_delt
    }
  }
  test_predicted <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt)
  train_predicted <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt)
  
  test_posterior <- pamr.predict(fit = training_model,newx = test_set$x,threshold = pam_delt,type = "posterior")
  train_posterior <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "posterior")
  
  ##correlation predicted
  #assign geneid and genenames for training model,
  #otherwise the pamr.listgenes would would throw error
  train_set$geneid <- 1:nrow(train_set$x) #make the $geneid components in train_set the same as training_model$gene.subset
  #if(!is.null(rownames(train_set$x))){
  train_set$genenames <- rownames(train_set$x) #make the $genenames components in train_set the same as training_model$gene.subset
  #}else{
  #  train_set$genename <- as.character(1:nrow(train_set$x))
  #}
  
  #pamr centroids
  pamr_centroids <- pamr.predict(fit = training_model,newx = train_set$x,threshold = pam_delt,type = "centroid")
  pam_info <- capture.output(pamr_survived_genes <- pamr.listgenes(fit = training_model,
                                                                   data = train_set,threshold = pam_delt,
                                                                   genenames = TRUE)[,"name"])
  
  #if you want to show the info from pamr.listgenes, uncomment the command below
  #cat(paste0(pam_info,collapse = "\n"))
  
  
  # pamr_survived_genes <- pamr.listgenes(fit = training_model,data = train_set,threshold = pam_delt,genenames = TRUE)[,"id"]
  if(length(pamr_survived_genes)<=1){
    message("Number of survived genes less or equal to 1 produced which would lead to error in the following analysis
            because structure of resulting data were changed!")
    return(list(delt=pam_delt,gene_numbers=length(pamr_survived_genes),traning_model=training_model,CR_model=cross_validation,
                expr=train_and_test_dataset))
  }
  pamr_centroids <- pamr_centroids[pamr_survived_genes,]
  
  #correlation
  pamr_centroids_test_cor <- cor(x = test_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
  pamr_centroids_train_cor <- cor(x = train_set$x[pamr_survived_genes,],y = pamr_centroids,method = "spearman")
  
  test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
  train_cor_predicted <- apply(pamr_centroids_train_cor,1,function(x)names(which.max(x)))
  
  predict_res<-list(train_set=data.frame(pam_predict=train_predicted,
                                         cor_predict=factor(train_cor_predicted,levels = levels(train_predicted)),
                                         stringsAsFactors = FALSE),#set both column of data.frame as factor, so that
                    #there would be the same columns in the table() function output, 
                    #otherwise there would be an error in caret::confusionMatrix
                    #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
                    #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
                    #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
                    test_set=data.frame(pam_predict=test_predicted,
                                        cor_predict=factor(test_cor_predicted,levels = levels(test_predicted)),
                                        stringsAsFactors = FALSE)#set both column of data.frame as factor, so that
                    #there would be the same columns in the table() function output, 
                    #otherwise there would be an error in caret::confusionMatrix
                    #try table(data.frame(a=c("a","b","c"),b=c("a","a","c"))) and 
                    #table(data.frame(a=factor(c("a","b","c"),levels=c("a","b","c")),
                    #                 b=factor(c("a","a","c"),levels=c("a","b","c"))))
  )
  kappas <- sapply(predict_res,function(x){
    confusionMatrix(table(x))$overall["Kappa"]
  })
  
  predict_res_error <- list(train_set=c(sum(train_predicted==train_and_test_dataset$train_set$y)/length(train_predicted),
                                        sum(train_cor_predicted==train_and_test_dataset$train_set$y)/length(train_cor_predicted)),
                            test_set=c(sum(test_predicted==train_and_test_dataset$test_set$y)/length(test_predicted),
                                       sum(test_cor_predicted==train_and_test_dataset$test_set$y)/length(test_cor_predicted)))
  predict_res_error$train_set <- setNames(predict_res_error$train_set,c("PAM_predicted","Cor_predicted"))
  predict_res_error$test_set <- setNames(predict_res_error$test_set,c("PAM_predicted","Cor_predicted"))
  
  if(!(all(train_and_test_dataset$train_set$y%in%levels(train_predicted))|
       all(train_and_test_dataset$test_set$y%in%levels(test_predicted))))
    stop("The levels of train or test predicted results is not match the real class")
  predict_res$train_set$real_class <- factor(train_and_test_dataset$train_set$y,levels = levels(train_predicted))
  predict_res$test_set$real_class <- factor(train_and_test_dataset$test_set$y,levels = levels(test_predicted))
  
  confusion_matrices <- list(train_set=list(PAM_and_real=confusionMatrix(table(predict_res$train_set[,c("pam_predict",
                                                                                                        "real_class")])),
                                            Cor_and_real=confusionMatrix(table(predict_res$train_set[,c("cor_predict",
                                                                                                        "real_class")]))),
                             test_set=list(PAM_and_real=confusionMatrix(table(predict_res$test_set[,c("pam_predict",
                                                                                                      "real_class")])),
                                           Cor_and_real=confusionMatrix(table(predict_res$test_set[,c("cor_predict",
                                                                                                      "real_class")]))))
  
  
  # list(kappas=kappas,predict_res=predict_res,train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor,
  #      centroid=pamr_centroids,delt=pam_delt,traning_and_testing=train_and_test_dataset)
  
  list(kappas=kappas,predict_res=predict_res,
       pamr_posterior=list(train_set=train_posterior,test_set=test_posterior),
       sample_cor=list(train_cor=pamr_centroids_train_cor,test_cor=pamr_centroids_test_cor),
       centroid=pamr_centroids,delt=pam_delt,traning_and_testing=train_and_test_dataset,
       error_rate=predict_res_error,predicted_confusions=confusion_matrices,
       CV_res=cross_validation,training_model=training_model)
}










#Updated function
#********************
#***nsccv_with_cor***
#********************

nsccv_with_cor <- function (x, y = NULL, proby = NULL, nfold = min(table(y)), 
                            folds = NULL, threshold = NULL, threshold.scale = NULL, 
                            survival.time = NULL, censoring.status = NULL, ngroup.survival = NULL, 
                            prior, object, ...) 
{
  this.call <- match.call()
  argy <- y
  if (is.null(y)) {
    y <- as.factor(apply(proby, 1, which.is.max))
  }
  n <- length(y)
  if (is.null(nfold) & is.null(survival.time)) {
    nfold <- min(table(y))
  }
  if (is.null(nfold) & !is.null(survival.time)) {
    nfold <- 10
  }
  if (is.null(survival.time)) {
    if (is.null(folds)) {
      folds <- pamr:::balanced.folds(y)
    }
  }
  if (!is.null(survival.time)) {
    if (is.null(folds)) {
      folds <- split(sample(1:n), rep(1:nfold, length = n))
    }
  }
  nfold <- length(folds)
  if (missing(prior)) {
    if (missing(object)) 
      prior <- table(y)/n
    else prior <- object$prior
  }
  if (missing(threshold)) {
    if (missing(object)) 
      stop("Must either supply threshold argument, or an nsc object")
    else {
      threshold <- object$threshold
      threshold.scale <- object$threshold.scale
      se.scale <- object$se.scale
    }
  }
  n.threshold <- length(threshold)
  yhat <- rep(list(y), n.threshold)
  names(yhat) <- paste(seq(n.threshold))
  yhat <- data.frame(yhat)
  
  # Cor_yhat <- matrix(0,nrow = ncol(x),ncol = n.threshold)
  # rownames(Cor_yhat) <- colnames(x)
  # colnames(Cor_yhat) <-  1:n.threshold
  Cor_yhat <- rep(list(y), n.threshold)
  names(Cor_yhat) <-  paste(seq(n.threshold))
  Cor_yhat <-  data.frame(Cor_yhat)
  
  
  n.class <- table(y)
  prob <- array(1, c(n, length(n.class), n.threshold))
  size <- double(n.threshold)
  hetero <- object$hetero
  cv.objects = vector("list", nfold)
  for (ii in 1:nfold) {
    cat("Fold", ii, ":")
    a <- pamr:::nsc(x[, -folds[[ii]]], y = argy[-folds[[ii]]], 
                    x[, folds[[ii]], drop = FALSE], proby = proby[-folds[[ii]], 
                    ], threshold = threshold, threshold.scale = threshold.scale, 
                    se.scale = se.scale, prior = prior, hetero = hetero, 
                    ..., remove.zeros = FALSE)
    size <- size + a$nonzero
    prob[folds[[ii]], , ] <- a$prob
    yhat[folds[[ii]], ] <- a$yhat
    cat("\n")
    cv.objects[[ii]] = a
    
    training_model <- pamr.train(data = list(x=x[, -folds[[ii]]],y=argy[-folds[[ii]]]),
                                 threshold = threshold)
    for(Nth_threshold in 1:length(threshold)){
      delt <- threshold[Nth_threshold]
      pamr_centroids <- pamr.predict(fit = training_model,newx = x[, folds[[ii]], drop = FALSE],threshold = delt,type = "centroid")
      
      pam_info <- capture.output(pamr_survived_genes <- pamr.listgenes(fit = training_model,
                                                                       data = list(x=x[, -folds[[ii]],drop=FALSE],y=argy[-folds[[ii]]],
                                                                                   geneid = 1:nrow(x),
                                                                                   genenames = rownames(x)),
                                                                       threshold = delt,
                                                                       genenames = TRUE)[,"name"])
      if(length(pamr_survived_genes)<=1){
        warning("Number of survived genes less or equal to 1, which lead to samples could not be predicted by Spearman Correlation!")
        Cor_yhat[folds[[ii]],Nth_threshold] <- NA
      }else{
        pamr_centroids <- pamr_centroids[pamr_survived_genes,]
        #correlation
        pamr_centroids_test_cor <- cor(x = x[pamr_survived_genes, folds[[ii]], drop = FALSE],y = pamr_centroids,method = "spearman")
        
        #test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x)names(which.max(x)))
        test_cor_predicted <- apply(pamr_centroids_test_cor,1,function(x){
          ifelse(all(is.na(x)),NA,names(which.max(x)))
        })
        
        
        
        # cat("Cor_yhat[folds[[ii]],Nth_threshold]: ",length(Cor_yhat[folds[[ii]],Nth_threshold]),"\n")
        # print(test_cor_predicted)
        # print(pamr_centroids_test_cor)
        # print(x[pamr_survived_genes, folds[[ii]], drop = FALSE])
        # print("")
        
        Cor_yhat[folds[[ii]],Nth_threshold] <- test_cor_predicted
      }
    }
    
  }
  if (missing(object)) 
    size <- round(size/nfold)
  else size <- object$nonzero
  error <- rep(NA, n.threshold)
  loglik <- error
  pvalue.survival <- error
  pvalue.survival.func <- function(group, survival.time, censoring.status, 
                                   ngroup.survival) {
    temp <- coxph(Surv(survival.time, censoring.status) ~ 
                    as.factor(group))
    loglik <- 2 * (temp$loglik[2] - temp$loglik[1])
    return(1 - pchisq(loglik, ngroup.survival - 1))
  }
  if (!is.null(proby)) {
    proby.temp <- proby
  }
  else if (!is.null(survival.time)) {
    proby.temp <- pamr.surv.to.class2(survival.time, censoring.status, 
                                      n.class = ngroup.survival)$prob
  }
  for (i in 1:n.threshold) {
    if (is.null(survival.time) & is.null(proby)) {
      error[i] <- sum(yhat[, i] != y)/n
    }
    if (!is.null(survival.time)) {
      temp <- c(yhat[, i], names(table(y)))
      Yhat <- model.matrix(~factor(temp) - 1, data = list(y = temp))
      Yhat <- Yhat[1:length(yhat[[ii]]), ]
      error[i] <- (length(yhat[, i]) - sum(Yhat * proby.temp))/n
    }
    if (is.null(survival.time)) {
      loglik[i] <- sum(log(prob[, , i][cbind(seq(1, n), 
                                             unclass(y))]))/n
    }
    if (!is.null(survival.time)) {
      pvalue.survival[i] <- pvalue.survival.func(yhat[, 
                                                      i], survival.time, censoring.status, ngroup.survival)
    }
  }
  obj <- list(threshold = threshold, error = error, loglik = loglik, 
              size = size, yhat = yhat, y = y, prob = prob, folds = folds, 
              cv.objects = cv.objects, pvalue.survival = pvalue.survival, 
              call = this.call,Cor_yhat=Cor_yhat)
  class(obj) <- "nsccv"
  obj
}




n_times_compare_train_test <- function(train_and_test_expr_with_label,n_times,train_proportion,auto_pam_delt,
                                       n_threshold,prior,...){
  compare_res <- list()
  i <- 1
  while(i <= n_times){
    message(paste0("The ",i,"th compare!"))
    #train_and_test_dat <- produce_train_test_set(expr_with_label = expr_with_label,train_proportion = train_proportion)
    # test_error <- PAM_test_error(train_and_test_dataset = train_and_test_dat,
    #                              auto_pam_delt = auto_pam_delt,n_threshold = n_threshold)
    train_and_test_dat <-  train_and_test_expr_with_label
    com_res <- compare_pam_and_correlation(train_and_test_dataset = train_and_test_dat,auto_pam_delt = auto_pam_delt,
                                           n_threshold = n_threshold,pamr_prior_class = prior,...)
    compare_res <- c(compare_res,list(com_res))
    i <- i+1
  }
  setNames(compare_res,paste0("The ",1:n_times,"st errors"))
}


#***************************************************
#***n_times_compare_train_test_for_Cor_validation***
#***************************************************

#This function is different from "n_times_compare_train_test"
#The output of this function have an additional element: output$`The nst errors`$CV_res$Cor_yhat
#The element "Cor_yhat" could not be generated by the function "compare_pam_and_correlation" used in "n_times_compare_train_test"
n_times_compare_train_test_for_Cor_validation <- function(train_and_test_expr_with_label,n_times,train_proportion,auto_pam_delt,
                                                          n_threshold,prior,...){
  compare_res <- list()
  i <- 1
  while(i <= n_times){
    message(paste0("The ",i,"th compare!"))
    #train_and_test_dat <- produce_train_test_set(expr_with_label = expr_with_label,train_proportion = train_proportion)
    # test_error <- PAM_test_error(train_and_test_dataset = train_and_test_dat,
    #                              auto_pam_delt = auto_pam_delt,n_threshold = n_threshold)
    train_and_test_dat <-  train_and_test_expr_with_label
    com_res <- compare_pam_and_correlation_for_cor(train_and_test_dataset = train_and_test_dat,auto_pam_delt = auto_pam_delt,
                                           n_threshold = n_threshold,pamr_prior_class = prior,...)
    compare_res <- c(compare_res,list(com_res))
    i <- i+1
  }
  setNames(compare_res,paste0("The ",1:n_times,"st errors"))
}

