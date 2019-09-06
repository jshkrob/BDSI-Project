
#######################------------------To run: Cmd F and replace with drug id label--------------------------#########################

# loading data
screening <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/screening.rds")
expr <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/expression.rds")
drugids <- c("1007", "1026", "1016", "1006", "1058", "1015", "1014", "1001", "1053",
             "1060", "1011", "1054", "1037", "1066", "1008")
library(readr)
library(pls)
library(e1071)
library(HiDimDA)
# packages:
library(tidyverse)
library(caret)
library(rpart)
library(randomForest)
library(naivebayes)
library(MASS)
library(glmnet)
library(caret)
library(rrlda)
library("class")

set.seed(23527)

colnames(screening)
View(screening)
drug=matrix(c(rep(0,times=236)),ncol=4)
d=unique(as.character(screening$DRUG_ID_lib))
d
library(naivebayes)
for(i in 1:59){
  x=subset(screening,DRUG_ID_lib==d[i])
  drug[i,1]=length(x[x$EFFECT=="TRUE",][,1])
  drug[i,2]=length(x[x$EFFECT=="FALSE",][,1])
  drug[i,3]=min((drug[i,1]/drug[i,2]),(drug[i,2]/drug[i,1]))
  drug[i,4]=as.numeric(d[i])
}
drug
View(drug)

findrug=drug[(drug[,3]>0.20)&((drug[,1]+drug[,2])>150),]
View(findrug)

#-------------------------------------------------------Functions:--------------------------------------------------------

# p-values 
cv_naive_bayes <- function(df, effect) {
  set.seed(23527)
  
  ### Naive Bayes CV w/ feature number:
  CV_NB <- data.frame(featurenum = numeric(10), error = numeric(10)) 
  for (i in 1:10) {
    n <- nrow(df)
    
    # training & test: 
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    df.train = df[train.ix,]
    effective.train = effect[train.ix]
    df.test = df[-train.ix,]
    effective.test = effect[-train.ix]
    
    best_q <- 0 
    CV_error <- vector(mode = "numeric")
    for(q in Q) {
      
      cv_error <- vector(mode = "numeric", length = 10)
      for(j in 1:10) {
        
        # splitting data into k folds:
        n2 <- nrow(df.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        df.train2 = df.train[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        df.validate = df.train[-train2.ix,]
        effective.validate = effective.train[-train2.ix]
        
        # Fit naive Bayes classifier on folds, compute CV error using validation set
        model_nb <- naive_bayes(x = df.train2[, 1:q], y = effective.train2)
        ytesthat <- predict(model_nb, df.validate[, 1:q])
        cv_error[j] <- 1 - mean(ytesthat == effective.validate)
        
      }
      # CV error for each feature number
      CV_error <- append(CV_error, mean(cv_error))
    }
    
    # Find q with the smallest error
    idx <- which.min(CV_error)
    best_q <- Q[idx]
    
    # Fit training data on selected feature number
    model_nb <- naive_bayes(x = df.train[, 1:best_q], y = effective.train)
    ytesthat <- predict(model_nb, df.test[, 1:q])
    test_error <- 1 - mean(ytesthat == effective.test)
    
    CV_NB[i, ] <- c(best_q, test_error)
  }
  # Out of 10 results, use the one w/ lowest cv error
  BEST_ERROR <- min(CV_NB$error)
  BEST_Q <- CV_NB$featurenum[which.min(CV_NB$error)]
  
  print(BEST_Q)
  return(c(BEST_ERROR, BEST_Q))
}
cv_lasso <- function(df, effect) {
  set.seed(23527)
  
  ### Lasso CV 
  CV_LASSO <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:2) {
      
      # FIRST, create TEST and TRAINING sets:
      n = nrow(df)
      train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
      df.train = df[train.ix,]
      effective.train = effect[train.ix]
      df.test = df[-train.ix,]
      effective.test = effect[-train.ix]
      
      # SECOND, select best # of features (q) given lambda via CV error
      best_q <- 0
      best_error <- 1
      
      for(q in Q) { 
        # CV step:
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
          # splitting data into k folds:
          n2 <- nrow(df.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          df.train2 = df.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          
          df.validate = df.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # train glm with regul. using q features w/ l_min choice & 19 training folds
          lasso.fit <- glmnet(x = as.matrix(df.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min)
          test_predic <- ifelse(predict(lasso.fit, as.matrix(df.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
          error_rate <- mean(test_predic != effective.validate)
          
          cv_error[j] <- error_rate
        }
        
        CV_error <- mean(cv_error)
        
        if (CV_error < best_error) {      
          best_q <- q
          best_error <- CV_error
        }
      }
      
      # THIRD, use feature number q to find best lamda for that q: 
      lasso.fit <- cv.glmnet(x = as.matrix(df.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 1)
      l_min <- lasso.fit$lambda.min
    }
    
    # use selected q and lamda on test set to estimtae error:
    lasso.fit <- glmnet(x = as.matrix(df.train[,1:best_q]), 
                        y = effective.train, 
                        family = "binomial", 
                        lambda = l_min)
    
    error <- mean((ifelse(predict(lasso.fit, as.matrix(df.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_LASSO[i,1] <- l_min
    CV_LASSO[i,2] <- best_q
    CV_LASSO[i,3] <- error
  }
  
  # return smallest error found among 10 examples:
  BEST_ERROR <- min(CV_LASSO$error)
  cat("error:", BEST_ERROR)
  cat("feature number:", best_q)
  cat("best parameters (lambda):", l_min)
  print(CV_LASSO)
  return(c(BEST_ERROR, best_q, l_min))
  
}
cv_ridge <- function(df, effect) {
  set.seed(23527)
  ### Ridge CV 
  CV_RIDGE <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:2) {
      
      # FIRST, create TEST and TRAINING sets:
      n = nrow(df)
      train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
      df.train = df[train.ix,]
      effective.train = effect[train.ix]
      df.test = df[-train.ix,]
      effective.test = effect[-train.ix]
      
      # SECOND, select best # of features (q) given lambda via CV error
      best_q <- 0
      best_error <- 1
      
      for(q in Q) { 
        # CV step:
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
          # splitting data into k folds:
          n2 <- nrow(df.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          df.train2 = df.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          
          df.validate = df.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # train glm with regul. using q features w/ l_min choice & 19 training folds
          lasso.fit <- glmnet(x = as.matrix(df.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min, alpha = 0)
          test_predic <- ifelse(predict(lasso.fit, as.matrix(df.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
          error_rate <- mean(test_predic != effective.validate)
          
          cv_error[j] <- error_rate
        }
        
        CV_error <- mean(cv_error)
        
        if (CV_error < best_error) {      
          best_q <- q
          best_error <- CV_error
        }
      }
      
      # THIRD, use feature number q to find best lamda for that q: 
      lasso.fit <- cv.glmnet(x = as.matrix(df.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 0)
      l_min <- lasso.fit$lambda.min
    }
    
    # use selected q and lamda on test set to estimtae error:
    lasso.fit <- glmnet(x = as.matrix(df.train[,1:best_q]), 
                        y = effective.train, 
                        family = "binomial", 
                        lambda = l_min)
    
    error <- mean((ifelse(predict(lasso.fit, as.matrix(df.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_RIDGE[i,1] <- l_min
    CV_RIDGE[i,2] <- best_q
    CV_RIDGE[i,3] <- error
  }
  
  # return smallest error found among 10 examples:
  BEST_ERROR <- min(CV_RIDGE$error)
  print(CV_RIDGE)
  cat("error:", BEST_ERROR)
  cat("feature number:", best_q)
  cat("best parameters (lambda):", l_min)
  return(c(BEST_ERROR, best_q, l_min))
  
}
cv_elas <- function(df, effect) {
  set.seed(23527)
  
  ### Elastic CV 
  CV_ELAS <- data.frame(alpha = numeric(10), lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  
  for(i in 1:10) {
    
    # FIRST, create TEST and TRAINING sets:
    n = nrow(df)
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    df.train = df[train.ix,]
    effective.train = effect[train.ix]
    df.test = df[-train.ix,]
    effective.test = effect[-train.ix]
    
    MCV <- vector(mode = "numeric")
    LAMBDA <- vector(mode = "numeric")
    
    # SECOND, select alpha and lambda based on CV:
    for(a in alpha) {
      # 10 fold cross validation for alpha = 0.0, 0.1, ..., 0.9, 1
      model_elas <- cv.glmnet(x = as.matrix(df.train), 
                              y = effective.train, family = "binomial",
                              type.measure = "class",
                              alpha = a)
      
      # smallest cv error and corresponding lambda
      min_lambda <- model_elas$lambda.min
      idx <- which(model_elas$lambda == min_lambda)
      cv_error <- model_elas$cvm[idx]
      
      # cv_error <- model_elas$cvm[which.min(model_elas$cvm)]
      # min_lambda <- model_elas$lambda[which.min(model_elas$cvm)]
      MCV <- append(MCV, cv_error)
      LAMBDA <- append(LAMBDA, min_lambda)
    }
    
    # find alpha that has lowest mcv error:
    idx = which.min(MCV)
    alpha_best <- alpha[idx]
    lambda_best <- LAMBDA[idx]
    
    # THRID, select best # of features (q) given lambda and alpha via CV error
    best_q <- 0
    best_error <- 1
    
    for(q in Q) { 
      
      # CV step:
      cv_error <- vector(mode = "numeric", length = 10)
      for(j in 1:10) {
        # splitting data into k folds:
        n2 <- nrow(df.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        df.train2 = df.train[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        df.validate = df.train[-train2.ix,]
        effective.validate = effective.train[-train2.ix]
        
        # Train using q features w/ l_min choice & 19 training folds
        elas.fit <- glmnet(x = as.matrix(df.train2[,1:q]), y = effective.train2, 
                           family = "binomial", lambda = lambda_best, alpha = alpha_best)
        
        test_predic <- ifelse(predict(elas.fit, as.matrix(df.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
        error_rate <- mean(test_predic != effective.validate)
        
        cv_error[j] <- error_rate
      }
      
      CV_error <- mean(cv_error)
      
      if (CV_error < best_error) {      
        best_q <- q
        best_error <- CV_error
      }
    }
    
    # Use selected q, alpha, and lamda on training set to estimate error:
    elas.fit <- glmnet(x = as.matrix(df.train[,1:best_q]), y = effective.train, 
                       family = "binomial", lambda = lambda_best, alpha = alpha_best)
    
    error <- mean((ifelse(predict(elas.fit, as.matrix(df.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_ELAS[i,1] <- alpha_best
    CV_ELAS[i,2] <- lambda_best
    CV_ELAS[i,3] <- best_q
    CV_ELAS[i,4] <- error
  }
  BEST_ERROR <- min(CV_ELAS$error)
  print(CV_ELAS)
  cat("error:", BEST_ERROR)
  cat("feature number:", best_q)
  cat("best parameters (alpha lambda):", alpha_best, lambda_best)
  return(c(BEST_ERROR, best_q, alpha_best, lambda_best))
  
}

# pca
cv_naive_bayes <- function(df, effect) {
  set.seed(23527)
  
  ### Naive Bayes CV w/ feature number:
  CV_NB <- data.frame(featurenum = numeric(10), error = numeric(10)) 
  for (i in 1:10) {
    n <- nrow(df)
    
    # training & test: 
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    df.train = df[train.ix,]
    effective.train = effect[train.ix]
    df.test = df[-train.ix,]
    effective.test = effect[-train.ix]
    
    best_q <- 0 
    CV_error <- vector(mode = "numeric")
    for(q in Q) {
      
      cv_error <- vector(mode = "numeric", length = 10)
      for(j in 1:10) {
        
        # splitting data into k folds:
        n2 <- nrow(df.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        df.train2 = df.train[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        df.validate = df.train[-train2.ix,]
        effective.validate = effective.train[-train2.ix]
        
        # Fit naive Bayes classifier on folds, compute CV error using validation set
        model_nb <- naive_bayes(x = df.train2[, 1:q], y = effective.train2)
        ytesthat <- predict(model_nb, df.validate[, 1:q])
        cv_error[j] <- 1 - mean(ytesthat == effective.validate)
        
      }
      # CV error for each feature number
      CV_error <- append(CV_error, mean(cv_error))
    }
    
    # Find q with the smallest error
    idx <- which.min(CV_error)
    best_q <- Q[idx]
    
    # Fit training data on selected feature number
    model_nb <- naive_bayes(x = df.train[, 1:best_q], y = effective.train)
    ytesthat <- predict(model_nb, df.test[, 1:q])
    test_error <- 1 - mean(ytesthat == effective.test)
    
    CV_NB[i, ] <- c(best_q, test_error)
  }
  # Out of 10 results, use the one w/ lowest cv error
  BEST_ERROR <- min(CV_NB$error)
  BEST_Q <- CV_NB$featurenum[which.min(CV_NB$error)]
  
  print(BEST_Q)
  return(c(BEST_ERROR, BEST_Q))
}
cv_lasso <- function(df, effect) {
  set.seed(23527)
  
  ### Lasso CV 
  CV_LASSO <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:2) {
      
      # FIRST, create TEST and TRAINING sets:
      n = nrow(df)
      train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
      df.train = df[train.ix,]
      effective.train = effect[train.ix]
      df.test = df[-train.ix,]
      effective.test = effect[-train.ix]
      
      # SECOND, select best # of features (q) given lambda via CV error
      best_q <- 0
      best_error <- 1
      
      for(q in Q) { 
        # CV step:
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
          # splitting data into k folds:
          n2 <- nrow(df.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          df.train2 = df.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          
          df.validate = df.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # train glm with regul. using q features w/ l_min choice & 19 training folds
          lasso.fit <- glmnet(x = as.matrix(df.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min)
          test_predic <- ifelse(predict(lasso.fit, as.matrix(df.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
          error_rate <- mean(test_predic != effective.validate)
          
          cv_error[j] <- error_rate
        }
        
        CV_error <- mean(cv_error)
        
        if (CV_error < best_error) {      
          best_q <- q
          best_error <- CV_error
        }
      }
      
      # THIRD, use feature number q to find best lamda for that q: 
      lasso.fit <- cv.glmnet(x = as.matrix(df.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 1)
      l_min <- lasso.fit$lambda.min
    }
    
    # use selected q and lamda on test set to estimtae error:
    lasso.fit <- glmnet(x = as.matrix(df.train[,1:best_q]), 
                        y = effective.train, 
                        family = "binomial", 
                        lambda = l_min)
    
    error <- mean((ifelse(predict(lasso.fit, as.matrix(df.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_LASSO[i,1] <- l_min
    CV_LASSO[i,2] <- best_q
    CV_LASSO[i,3] <- error
  }
  
  # return smallest error found among 10 examples:
  BEST_ERROR <- min(CV_LASSO$error)
  cat("error:", BEST_ERROR)
  cat("feature number:", best_q)
  cat("best parameters (lambda):", l_min)
  print(CV_LASSO)
  return(c(BEST_ERROR, best_q, l_min))
  
}
cv_ridge <- function(df, effect) {
  set.seed(23527)
  ### Ridge CV 
  CV_RIDGE <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:2) {
      
      # FIRST, create TEST and TRAINING sets:
      n = nrow(df)
      train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
      df.train = df[train.ix,]
      effective.train = effect[train.ix]
      df.test = df[-train.ix,]
      effective.test = effect[-train.ix]
      
      # SECOND, select best # of features (q) given lambda via CV error
      best_q <- 0
      best_error <- 1
      
      for(q in Q) { 
        # CV step:
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
          # splitting data into k folds:
          n2 <- nrow(df.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          df.train2 = df.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          
          df.validate = df.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # train glm with regul. using q features w/ l_min choice & 19 training folds
          lasso.fit <- glmnet(x = as.matrix(df.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min, alpha = 0)
          test_predic <- ifelse(predict(lasso.fit, as.matrix(df.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
          error_rate <- mean(test_predic != effective.validate)
          
          cv_error[j] <- error_rate
        }
        
        CV_error <- mean(cv_error)
        
        if (CV_error < best_error) {      
          best_q <- q
          best_error <- CV_error
        }
      }
      
      # THIRD, use feature number q to find best lamda for that q: 
      lasso.fit <- cv.glmnet(x = as.matrix(df.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 0)
      l_min <- lasso.fit$lambda.min
    }
    
    # use selected q and lamda on test set to estimtae error:
    lasso.fit <- glmnet(x = as.matrix(df.train[,1:best_q]), 
                        y = effective.train, 
                        family = "binomial", 
                        lambda = l_min)
    
    error <- mean((ifelse(predict(lasso.fit, as.matrix(df.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_RIDGE[i,1] <- l_min
    CV_RIDGE[i,2] <- best_q
    CV_RIDGE[i,3] <- error
  }
  
  # return smallest error found among 10 examples:
  BEST_ERROR <- min(CV_RIDGE$error)
  print(CV_RIDGE)
  cat("error:", BEST_ERROR)
  cat("feature number:", best_q)
  cat("best parameters (lambda):", l_min)
  return(c(BEST_ERROR, best_q, l_min))
  
}
cv_elas <- function(df, effect) {
  set.seed(23527)
  
  ### Elastic CV 
  CV_ELAS <- data.frame(alpha = numeric(10), lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  
  for(i in 1:10) {
    
    # FIRST, create TEST and TRAINING sets:
    n = nrow(df)
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    df.train = df[train.ix,]
    effective.train = effect[train.ix]
    df.test = df[-train.ix,]
    effective.test = effect[-train.ix]
    
    MCV <- vector(mode = "numeric")
    LAMBDA <- vector(mode = "numeric")
    
    # SECOND, select alpha and lambda based on CV:
    for(a in alpha) {
      # 10 fold cross validation for alpha = 0.0, 0.1, ..., 0.9, 1
      model_elas <- cv.glmnet(x = as.matrix(df.train), 
                              y = effective.train, family = "binomial",
                              type.measure = "class",
                              alpha = a)
      
      # smallest cv error and corresponding lambda
      min_lambda <- model_elas$lambda.min
      idx <- which(model_elas$lambda == min_lambda)
      cv_error <- model_elas$cvm[idx]
      
      # cv_error <- model_elas$cvm[which.min(model_elas$cvm)]
      # min_lambda <- model_elas$lambda[which.min(model_elas$cvm)]
      MCV <- append(MCV, cv_error)
      LAMBDA <- append(LAMBDA, min_lambda)
    }
    
    # find alpha that has lowest mcv error:
    idx = which.min(MCV)
    alpha_best <- alpha[idx]
    lambda_best <- LAMBDA[idx]
    
    # THRID, select best # of features (q) given lambda and alpha via CV error
    best_q <- 0
    best_error <- 1
    
    for(q in Q) { 
      
      # CV step:
      cv_error <- vector(mode = "numeric", length = 10)
      for(j in 1:10) {
        # splitting data into k folds:
        n2 <- nrow(df.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        df.train2 = df.train[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        df.validate = df.train[-train2.ix,]
        effective.validate = effective.train[-train2.ix]
        
        # Train using q features w/ l_min choice & 19 training folds
        elas.fit <- glmnet(x = as.matrix(df.train2[,1:q]), y = effective.train2, 
                           family = "binomial", lambda = lambda_best, alpha = alpha_best)
        
        test_predic <- ifelse(predict(elas.fit, as.matrix(df.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
        error_rate <- mean(test_predic != effective.validate)
        
        cv_error[j] <- error_rate
      }
      
      CV_error <- mean(cv_error)
      
      if (CV_error < best_error) {      
        best_q <- q
        best_error <- CV_error
      }
    }
    
    # Use selected q, alpha, and lamda on training set to estimate error:
    elas.fit <- glmnet(x = as.matrix(df.train[,1:best_q]), y = effective.train, 
                       family = "binomial", lambda = lambda_best, alpha = alpha_best)
    
    error <- mean((ifelse(predict(elas.fit, as.matrix(df.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_ELAS[i,1] <- alpha_best
    CV_ELAS[i,2] <- lambda_best
    CV_ELAS[i,3] <- best_q
    CV_ELAS[i,4] <- error
  }
  BEST_ERROR <- min(CV_ELAS$error)
  print(CV_ELAS)
  cat("error:", BEST_ERROR)
  cat("feature number:", best_q)
  cat("best parameters (alpha lambda):", alpha_best, lambda_best)
  return(c(BEST_ERROR, best_q, alpha_best, lambda_best))
  
}

#------------------------------------------------------DRUG ID: 1007------------------------------------------------------------------

effective.1007 <- subset(screening, DRUG_ID_lib == "1007")[,c("CL","EFFECT")]
expr.1007 <- expr[as.character(effective.1007$CL),]
Y=effective.1007$EFFECT
Y
dim(effective.1007)
dim(expr.1007)

get.p <- function(dat, labs){
  # split the data into effective and ineffective
  effect <- dat[labs]
  ineffect <- dat[!labs]
  
  # calculate the two sample means
  effect.bar <- mean(effect)
  ineffect.bar <- mean(ineffect)
  
  # calculate the two sample variances
  v.effect <- var(effect)
  v.ineffect <- var(ineffect)
  
  # calculate the sample sizes
  n.effect <- length(effect)
  n.ineffect <- length(ineffect)
  
  # calculate the sd
  s_pooled = (n.effect*v.effect +n.ineffect*v.ineffect)/(n.effect+n.ineffect-2)
  s <- sqrt(s_pooled*(1/n.effect + 1/n.ineffect))
  
  # calculate the test statistic
  T_abs <- abs((effect.bar - ineffect.bar)/s)
  pval = 2*(1-pt(T_abs,df = n.effect+n.ineffect - 2))
  return(pval)
}

############################################################ Naive Bayes

n=nrow(expr.1007)
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)

# t-test feature selection: # of features from 50 to 500
dim_seq=seq(50,500,50)
dim_seq

# 5 trials:
for(i in 1:5){
  
  # Select training and testing sets: 4-1 split
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = expr.1007[train.ix,]
  Ytrain = effective.1007$EFFECT[train.ix]
  Xtest = expr.1007[-train.ix,]
  Ytest = effective.1007$EFFECT[-train.ix]
  
  # Average CV error for each m in dim_seq
  dim.error <- vector(mode = "numeric")
  
  # find CV error for each feature number
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    # Average missiclassification error over 5 trials
    for(j in 1:5){
      
      # Create 5 folds: 
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],m)
      
      #cv.features[j,] = pvals_sel$GENE
      
      # selecting the top *m* features (from dim_seq)
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      
      model1 = naive_bayes(x=Xtrain2,y =Ytrain2)
      Ypred = predict(model1,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  
  # feature number that produced minimum CV error:
  idx <- which.min(dim.error)
  best.dim[i] <- dim_seq[idx]
  
  # important genes from training set: using best feature #...
  pvals <- data.frame(GENE = colnames(expr.1007))
  pvals$p = apply(Xtrain,2,get.p, Ytrain)
  pvals_sel = head(pvals[order(pvals$p),],best.dim[i])
  Xtrain1 = Xtrain[,pvals_sel$GENE]
  Xtest1 = Xtest[,pvals_sel$GENE]
  
  model2= naive_bayes(x=Xtrain1,y = Ytrain)
  Ypred = predict(model2,Xtest1)
  cv.error[i]=mean(Ypred!= Ytest)
}

# Test errors after Naive Bayes: best feature #, CV for feature #
tab_naivebayes= cbind(best.dim,cv.error)[1:5,]
tab_naivebayes
testerror=c(rep(0,times=10))
testerror[1]=(cv.error[1]+cv.error[4])/2
testerror

############################################################# LDA:
n=nrow(expr.1007)
##cv.error <- vector(mode = "numeric", length = 5)
##best.dim <- vector(mode = "numeric", length = 5)
dim_seq=seq(50,500,50)
dim_seq

for(i in 1:5){
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = expr.1007[train.ix,]
  Ytrain = effective.1007$EFFECT[train.ix]
  Xtest = expr.1007[-train.ix,]
  Ytest = effective.1007$EFFECT[-train.ix]
  
  dim.error <- vector(mode = "numeric")
  
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      # t-test feature selection:
      
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],m)
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      
      model1 = lda(Xtrain2,Ytrain2)
      Ypred = predict(model1,Xvalidate)$class
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim[i] <- dim_seq[idx]
  pvals <- data.frame(GENE = colnames(expr.1007))
  pvals$p = apply(Xtrain,2,get.p, Ytrain)
  pvals_sel = head(pvals[order(pvals$p),],best.dim[i])
  Xtrain1 = Xtrain[,pvals_sel$GENE]
  Xtest1 = Xtest[,pvals_sel$GENE]
  
  model2= lda(Xtrain1,Ytrain)
  Ypred = predict(model2,Xtest1)$class
  cv.error[i]=mean(Ypred!= Ytest)
}
tab_lda= cbind(best.dim,cv.error)[1:5,]
tab_lda
testerror[2]=(cv.error[3]+cv.error[4]+cv.error[5])/3
testerror

################################################## GLMNET
library(glmnet)

dim_seq=seq(100,500,50)

# Authentic lambdas after 1 model fit:
lambda = cv.glmnet(as.matrix(expr.1007),as.double(effective.1007$EFFECT))$lambda
lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
n <- nrow(expr.1007)

# create a vector to hold the 10 error rate

# create a vector to hold the 10 best values of k
##best.lambda <-matrix(0,nrow=9,ncol=10)
##finerror= vector(mode = "numeric", length = 5)
##finlambda= vector(mode = "numeric", length = 5)
##findim= vector(mode = "numeric", length = 5)
#best.features <- matrix(NA,nrow = 10,ncol = 100)

###for(i in 1:5){
train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]

# MCV Matrix (feature # v. lambda)
lambda.error <- matrix(0,nrow=9,ncol=10)
lambda.error

# CV for each feature number and lambda:
for(k in 1:9){
  
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],dim_seq[k])
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      
      # fit glm w/ certain lambda and feature number
      fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = lambda_seq[m])
      Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error[k,m] = mean(inside.error)
  }
}

# find the k with the smallest error

#-------------------------------  I think it should be idx <- which.min(lambda.error)
#idx <- which.min(cost.error)
idx <- which.min(lambda.error)
if(idx%%9!=0) {
  best_dim=dim_seq[idx%%9]
  best_cost=cost_seq[(idx%/%9+1)]
} else {
  best_dim=dim_seq[9]
  best_cost=cost_seq[(idx%/%9)]
}

# Using best_dim and best_cost: 
##best.lambda[i,k] <- lambda_seq[idx]
pvals <- data.frame(GENE = colnames(expr.1007))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]
  
fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda )
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_glmnet=c(best_lambda,best_dim,cv.error)
tab_glmnet  

########################################################### KNN
library(caret)

# create the vector of tuning parameters
dim_seq=seq(100,500,50)
lambda_seq=seq(1,19,2)
n <- nrow(expr.1007)

train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)
  
lambda.error <- matrix(0,nrow=9,ncol=10)
lambda.error
for(i in 1:9){
  for(m in 1:10){
      
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
        
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],dim_seq[i])
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
        
      model5a = knn3(Ytrain2 ~ ., data=train2,  k = lambda_seq[m]) 
      Ypred = ifelse(predict(model5a,Xvalidate)[,2]>=.5,1,0)
      inside.error[j] <- mean(Ypred != Yvalidate)
        
    }
      
    lambda.error[i,m] = mean(inside.error)
  }
}
# find the k with the smallest error
#---------------------------------------------------------- same error?
#idx <- which.min(cost.error)
idx <- which.min(lambda.error)
if(idx%%9!=0){
  best_dim=dim_seq[idx%%9]
  best_lambda=lambda_seq[(idx%/%9+1)]
} else {
  best_dim=dim_seq[9]
  best_lambda=lambda_seq[(idx%/%9)]
}
pvals <- data.frame(GENE = colnames(expr.1007))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)
  
model5a = knn3(Ytrain ~ ., data=train1,  k = best_lambda) 
Ypred = ifelse(predict(model5a,Xtest1)[,2]>=.5,1,0)
cv.error <- mean(Ypred!= Ytest)
tab_knn=c(best_lambda,best_dim,cv.error)
tab_knn
warnings()
#######################RandomForest

library(randomForest)
  
n=nrow(expr.1007)
dim_seq=seq(500,2000,500)
dim_seq
  
    
train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
dim.error <- vector(mode = "numeric")
    
for(m in 1:4) {
      
  inside.error <- vector(mode = "numeric", length = 5)    
  #cv.features <- matrix(NA,nrow = 5,ncol = 100)
      
  for(j in 1:5) {
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
        
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],dim_seq[m])
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
        
      model6 <- randomForest(as.factor(Ytrain2) ~ ., data=train2)
      Ypred = predict(model6,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
  }
    
  dim.error[m] <- mean(inside.error)
}

# Find the minimum error in dim.error

idx <- which.min(dim.error)
best.dim<- dim_seq[idx]
pvals <- data.frame(GENE = colnames(expr.1007))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best.dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)
model6 <- randomForest(as.factor(Ytrain) ~ ., data=train1)
Ypred = predict(model6,Xtest1)
cv.error=mean(Ypred!= Ytest)
cv.error
tab_rf=c(best.dim,cv.error)
tab_rf

############################################# SVM
library(e1071)

dim_seq=seq(50,500,50)
cost_seq=seq(0.25,2.5,0.25)
cost_seq
# create the vector of tuning parameters
n <- nrow(expr.1007)

train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]

cost.error <- matrix(0,nrow=10,ncol=10)
cost.error
for(i in 1:10){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],dim_seq[i])
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      
      
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="linear",cost=cost_seq[m])
      Ypred<-predict(svm_model,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    cost.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(cost.error)
if(idx%%10!=0) { 
  best_dim=dim_seq[idx%%10]
  best_cost=cost_seq[(idx%/%10+1)]
} else {
  best_dim=dim_seq[10]
  best_cost=cost_seq[(idx%/%10)]
}
pvals <- data.frame(GENE = colnames(expr.1007))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]

# Compute test error for SVM
svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="linear",cost=best_cost)
Ypred<-predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svm=c(best_cost,best_dim,cv.error)
tab_svm

###################################################### GLMNET-RIDGE

library(glmnet)

# Fit to find basic lambda values
dim_seq=seq(100,500,50)
lambda = cv.glmnet(as.matrix(expr.1007),as.double(effective.1007$EFFECT),alpha=0)$lambda
lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
n <- nrow(expr.1007)
lambda


# create a vector to hold the 10 error rate

# create a vector to hold the 10 best values of k

train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]

lambda.error <- matrix(0,nrow=9,ncol=10)
lambda.error
for(k in 1:9){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],dim_seq[k])
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      
      fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = lambda_seq[m],alpha=0)
      Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error[k,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(cost.error)
if(idx%%9!=0) {
  best_dim=dim_seq[idx%%9]
  best_cost=cost_seq[(idx%/%9+1)]
} else {
  best_dim=dim_seq[9]
  best_cost=cost_seq[(idx%/%9)]
}
##best.lambda[i,k] <- lambda_seq[idx]
pvals <- data.frame(GENE = colnames(expr.1007))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda ,alpha=0)
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_glmnetridge=c(best_lambda,best_dim,cv.error)
tab_glmnetridge

#################################################### SVM-POLYNOMIAL 
  
library(e1071)

# create the vector of tuning parameters
dim_seq=seq(50,500,50)
cost_seq=seq(0.25,2.5,0.25)

n <- nrow(expr.1007)

train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]

cost.error <- matrix(0,nrow=10,ncol=10)
cost.error

for(i in 1:10){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      pvals <- data.frame(GENE = colnames(expr.1007))
      pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
      pvals_sel = head(pvals[order(pvals$p),],dim_seq[i])
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,pvals_sel$GENE]
      Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
      Yvalidate = Ytrain[-train2.ix]
      
      
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="polynomial",cost=cost_seq[m])
      Ypred<-predict(svm_model,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    cost.error[i,m] =mean(inside.error)
  }
}

# find the k with the smallest error
idx <- which.min(cost.error)
if(idx%%10!=0) {
  best_dim=dim_seq[idx%%10]
  best_cost=cost_seq[(idx%/%10+1)]
} else {
  best_dim=dim_seq[10]
  best_cost=cost_seq[(idx%/%10)]
}

pvals <- data.frame(GENE = colnames(expr.1007))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]

svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="polynomial",cost=best_cost)
Ypred<-predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svmpoly=c(best_cost,best_dim,cv.error)
tab_svmpoly

################################################################## RDA
  
n=nrow(expr.1007)
library(rda)

train.ix <- sample(1:n, 0.6*n)
Xtrain = expr.1007[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = expr.1007[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train=data.frame(Xtrain,Ytrain)

fit<-rda(x = t(as.matrix(Xtrain)),y = as.double(Ytrain), xnew = t(Xtest),ynew = as.double(Ytest),delta = 0.0001,alpha = 0)

# Ypred<-predict(fit,t(Xtrain),Ytrain,t(Xtest))
# error<-mean(Ypred!=Ytest)
# error

#------------------------------------------------------Feature Selection: PCA------------------------------------------------------------------

#######################################################  PCA

#expr=readRDS("expression.rds")
#screening = readRDS("screening.rds")
effective.1007 <- subset(screening, DRUG_ID_lib == "1001")[,c("CL","CELL_LINE_NAME","EFFECT")]
dim(effective.1007)
expr.1007 <- expr[as.character(effective.1007$CL),]
center.expr=scale(expr.1007,scale=F)
View(center.expr)
pca <- svd(center.expr, nu = 20, nv = 0)  
pca  
m = sum(cumsum(pca$d^2)/sum(pca$d^2)<0.9)
m

# Extraction of all pc components:
pca_comp=prcomp(center.expr)
pc.scoresdat<-pca_comp$x
dim(pc.scoresdat)
pc.scores <- cbind(effective.1007, pc.scoresdat)
dim(pc.scores)
saveRDS(pc.scores,"pcaexpression-90%variance.rds")

######################################################## PCR

library(pls)

n=dim(center.expr)[1]
n
dim_seq=seq(5,40,5)
dim_seq
get.p <- function(dat, labs){
  # split the data into effective and ineffective
  effect <- dat[labs]
  ineffect <- dat[!labs]
  
  # calculate the two sample means
  effect.bar <- mean(effect)
  ineffect.bar <- mean(ineffect)
  
  # calculate the two sample variances
  v.effect <- var(effect)
  v.ineffect <- var(ineffect)
  
  # calculate the sample sizes
  n.effect <- length(effect)
  n.ineffect <- length(ineffect)
  
  # calculate the sd
  s_pooled = (n.effect*v.effect +n.ineffect*v.ineffect)/(n.effect+n.ineffect-2)
  s <- sqrt(s_pooled*(1/n.effect + 1/n.ineffect))
  
  # calculate the test statistic
  T_abs <- abs((effect.bar - ineffect.bar)/s)
  pval = 2*(1-pt(T_abs,df = n.effect+n.ineffect - 2))
  return(pval)
}
  
pvals <- data.frame(GENE = colnames(center.expr))
pvals$p = apply(center.expr,2,get.p, Y)
pvals_sel = head(pvals[order(pvals$p),],2000)
center.expr1000=center.expr[,pvals_sel$GENE]
  
train.ix <- sample(1:n, (4/5)*n)
Xtrain = center.expr1000[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = center.expr1000[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train.dat <- data.frame(class = Ytrain, Xtrain)
test.dat <- data.frame(class =Ytest, Xtest)
dim.error <- vector(mode = "numeric")
  
# CV for number of componenets in PCR
for(m in dim_seq){
    
  inside.error <- vector(mode = "numeric", length = 5)
  #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
  for(j in 1:5){
    n2 <- nrow(Xtrain)
    train2.ix <- sample(1:n2, (4/5)*n2)
    Xtrain2 = Xtrain[train2.ix,]
    Ytrain2 = Ytrain[train2.ix]
    train2.dat <- data.frame(class = Ytrain2, Xtrain2)
      
    Xvalidate = Xtrain[-train2.ix,]
    Yvalidate = Ytrain[-train2.ix]
    test2.dat <- data.frame(class =Yvalidate, Xvalidate)
      
      
    pcr.fit <- pcr(class ~ .,data=train2.dat,scale = TRUE, ncomp = m)
      
    # now use that fit to get predictions for the testing data
    pcr.pred <- predict(pcr.fit,test2.dat,ncomp=m)
    Ypred <- ifelse(pcr.pred >= 0.5,1,0)
    inside.error[j] <- mean(Ypred != Yvalidate)
  }
    
  dim.error <-append(dim.error,mean(inside.error))
}

idx <- which.min(dim.error)
best.dim<- dim_seq[idx]
pcr.fit <- pcr(class~ .,data=train.dat,scale = TRUE, ncomp =best.dim)
  
# now use that fit to get predictions for the testing data
pcr.pred <- predict(pcr.fit,test.dat,ncomp=best.dim)
Ypred <- ifelse(pcr.pred >= 0.5,1,0)
cv.error=mean(Ypred!= Ytest)

# Test error for PCR model:
tab_pcr=cbind(best.dim,cv.error)

tab_pcr

#################################################### PC-LDA

n=nrow(expr.1007)
pca_comp=prcomp(center.expr)
pc.scoresdat<-pca_comp$x
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)

dim_seq=seq(5,40,5)
library(MASS)
for(i in 1:5)
{train.ix <- sample(1:n, 0.8*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]

dim.error <- vector(mode = "numeric")
for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,1:m]
      Ytrain2 = Ytrain[train2.ix]
      Xvalidate = Xtrain[-train2.ix,1:m]
      Yvalidate = Ytrain[-train2.ix]
      
      model1 = lda(Xtrain2,Ytrain2)
      Ypred = predict(model1,Xvalidate)$class
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim[i]<- dim_seq[idx]
 
  Xtrain1 = Xtrain[,1:best.dim[i]]
  Xtest1 = Xtest[,1:best.dim[i]]
  
  model2= lda(Xtrain1,Ytrain)
  Ypred = predict(model2,Xtest1)$class
  cv.error[i]=mean(Ypred!= Ytest)
}
tab_pclda= cbind(best.dim,cv.error)[1:5,]
tab_pclda


################################################### PCA-KNN
library(caret)

dim_seq=seq(5,40,5)
lambda_seq=seq(1,19,2)
# create the vector of tuning parameters
n <- nrow(expr.1007)
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

lambda.error <- matrix(0,nrow=8,ncol=10)
lambda.error
for(i in 1:8){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,1:dim_seq[i]]
      Xvalidate = Xtrain[-train2.ix,1:dim_seq[i]]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      
      model5a = knn3(Ytrain2 ~ ., data=train2,  k = lambda_seq[m]) 
      Ypred = ifelse(predict(model5a,data.frame(Xvalidate))[,2]>=.5,1,0)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(lambda.error)
if(idx%%8!=0)
{best_dim=dim_seq[idx%%8]
best_lambda=lambda_seq[(idx%/%8+1)]
}else
{best_dim=dim_seq[8]
best_lambda=lambda_seq[(idx%/%8)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

# Compute test error:
model5a = knn3(Ytrain ~ ., data=train1,  k = best_lambda) 
Ypred = ifelse(predict(model5a,data.frame(Xtest1))[,2]>=.5,1,0)
cv.error <- mean(Ypred!= Ytest)
tab_knn=c(best_lambda,best_dim,cv.error)
tab_knn


################################################### PCA-SVM Linear
library(e1071)

dim_seq=seq(5,40,5)
cost_seq=seq(0.25,2.5,0.25)
# create the vector of tuning parameters
n <- nrow(pc.scoresdat)
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

cost.error <- matrix(0,nrow=8,ncol=10)
cost.error
for(i in 1:8){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,1:dim_seq[i]]
      Xvalidate = Xtrain[-train2.ix,1:dim_seq[i]]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="linear",cost=cost_seq[m])
      Ypred = predict(svm_model,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    cost.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(cost.error)
if(idx%%8!=0)
{best_dim=dim_seq[idx%%8]
best_cost=cost_seq[(idx%/%8+1)]
}else
{best_dim=dim_seq[8]
best_cost=cost_seq[(idx%/%8)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="linear",cost=best_cost)
Ypred = predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svm=c(best_cost,best_dim,cv.error)
tab_svm


################################################### PCA-SVM Radial
library(e1071)

dim_seq=seq(5,40,5)
cost_seq=seq(0.25,2.5,0.25)
# create the vector of tuning parameters
n <- nrow(pc.scoresdat)
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

cost.error <- matrix(0,nrow=8,ncol=10)
cost.error
for(i in 1:8){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,1:dim_seq[i]]
      Xvalidate = Xtrain[-train2.ix,1:dim_seq[i]]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="radial",cost=cost_seq[m])
      Ypred = predict(svm_model,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    cost.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(cost.error)
if(idx%%8!=0)
{best_dim=dim_seq[idx%%8]
best_cost=cost_seq[(idx%/%8+1)]
}else
{best_dim=dim_seq[8]
best_cost=cost_seq[(idx%/%8)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="radial",cost=best_cost)
Ypred = predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svm=c(best_cost,best_dim,cv.error)
tab_svm


################################################### PCA-SVM Polynomial
library(e1071)

dim_seq=seq(5,40,5)
cost_seq=seq(0.25,2.5,0.25)
# create the vector of tuning parameters
n <- nrow(pc.scoresdat)
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

cost.error <- matrix(0,nrow=8,ncol=10)
cost.error
for(i in 1:8){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,1:dim_seq[i]]
      Xvalidate = Xtrain[-train2.ix,1:dim_seq[i]]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="polynomial",cost=cost_seq[m])
      Ypred = predict(svm_model,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    cost.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(cost.error)
if(idx%%8!=0)
{best_dim=dim_seq[idx%%8]
best_cost=cost_seq[(idx%/%8+1)]
}else
{best_dim=dim_seq[8]
best_cost=cost_seq[(idx%/%8)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="polynomial",cost=best_cost)
Ypred = predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svm=c(best_cost,best_dim,cv.error)
tab_svm

####################################################### PCA-Naive Bayes

n=nrow(expr.1007)
pca_comp=prcomp(center.expr)
pc.scoresdat<-pca_comp$x
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)

dim_seq=seq(5,40,5)
library(naivebayes)
for(i in 1:5)
{train.ix <- sample(1:n, 0.8*n)
  Xtrain = pc.scoresdat[train.ix,]
  Ytrain = effective.1007$EFFECT[train.ix]
  Xtest = pc.scoresdat[-train.ix,]
  Ytest = effective.1007$EFFECT[-train.ix]
  
  dim.error <- vector(mode = "numeric")
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,1:m]
      Ytrain2 = Ytrain[train2.ix]
      Xvalidate = Xtrain[-train2.ix,1:m]
      Yvalidate = Ytrain[-train2.ix]
      
      model1 = naive_bayes(x=Xtrain2,y =Ytrain2)
      Ypred = predict(model1,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim[i]<- dim_seq[idx]
  
  Xtrain1 = Xtrain[,1:best.dim[i]]
  Xtest1 = Xtest[,1:best.dim[i]]
  
  model2= naive_bayes(x=Xtrain1,y = Ytrain)
  Ypred = predict(model2,Xtest1)
  cv.error[i]=mean(Ypred!= Ytest)
}
tab_pcnb= cbind(best.dim,cv.error)[1:5,]
tab_pcnb


################################################### PCA-Ridge
library(glmnet)

dim_seq=seq(5,40,5)
lambda = cv.glmnet(as.matrix(pc.scoresdat),as.double(effective.1007$EFFECT),alpha=0)$lambda

lambda_seq = seq(min(lambda),max(lambda),length.out = 10) 
lambda_seq
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1007$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1007$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

lambda.error <- matrix(0,nrow=8,ncol=10)
lambda.error
for(i in 1:8){
  for(m in 1:10){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      
      #cv.features[j,] = pvals_sel$GENE
      Xtrain2 = Xtrain2[,1:dim_seq[i]]
      Xvalidate = Xtrain[-train2.ix,1:dim_seq[i]]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      
      fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = lambda_seq[m],alpha=0)
      Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(lambda.error)
if(idx%%8!=0)
{best_dim=dim_seq[idx%%8]
best_lambda=lambda_seq[(idx%/%8+1)]
}else
{best_dim=dim_seq[8]
best_lambda=lambda_seq[(idx%/%8)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda ,alpha=0)
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_ridge=c(best_lambda,best_dim,cv.error)
tab_ridge

################################################ PCA Ridge, Lasso, Elas w/ Iterations
Q <- seq(10, col_number(pc.scoresdat), 5)

cv_ridge(as.data.frame(pc.scoresdat), effective.1007$EFFECT)
cv_lasso(as.data.frame(pc.scoresdat), effective.1007$EFFECT)
cv_elas(as.data.frame(pc.scoresdat), effective.1007$EFFECT)

############################################################################ New algorithm:

# change of basis into the proper principal component
test %*% pca$rotation






































