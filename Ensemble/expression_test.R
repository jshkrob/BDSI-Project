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

# modeling functions

#Same t-test function as we used for gene expression data
#Params: dat is the dataframe, labs are the labels 
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
  #s <- sqrt((v.effect/n.effect) + (v.ineffect/n.ineffect))
  
  # calculate the test statistic
  T_abs <- abs((effect.bar - ineffect.bar)/s)
  pval = 2*(1-pt(T_abs,df = n.effect+n.ineffect - 2))
  return(pval)
}

#Params: df, the data 
#effect: the classifiers for whether effective or ineffective
rand_forest <- function(df,effect) {# split the centered data into training and testing
  set.seed(23527)
  errors <- vector(mode='numeric', length=10)
  
  n <- dim(df)[1]
  
  
  # training set
  for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- df[train.index,]
    y.train <- effect[train.index]
    train.dat <- data.frame(class = y.train, x.train)
    dim(train.dat)
    # testing set
    x.test <- df[-train.index,]
    y.test <- effect[-train.index]
    test.dat <- data.frame(class = y.test, x.test)
    dim(test.dat)
    rf.fit <- randomForest(as.factor(class) ~ ., data=train.dat)
    #Ytrainhat = rf.fit$predicted
    preds <- predict(rf.fit,x.test)
    #trainingerror$rf = 1 - mean(Ytrainhat==Ytrain)
    #testerror$rf = 1 - mean(Ytesthat==Ytest)
    
    # get the classification success rate
    errors[i] <- mean(preds != y.test)    
  }
  
  return(mean(errors))
}

#New M param helps specify the sequences needed for your param
#sequence step must be 1 in this case but i'll fix it later
cv.rand_forest <- function(df,effect, M) {
  
  
  n <- nrow(df)
  error <- vector(mode='numeric', length=length(M))
  
  train.index <- sample(1:n, round(.9*n))
  x.train <- df[train.index,]
  y.train <- effect[train.index]
  train.dat <- data.frame(class = y.train, x.train)
  dim(train.dat)
  
  # testing set
  x.test <- df[-train.index,]
  y.test <- effect[-train.index]
  test.dat <- data.frame(class = y.test, x.test)
  
  for (i in M) {
    error[((i-M[1])/(M[2]-M[1]))+1] <- rand_forest(x.train[,1:i],y.train)
  }
  
  cbind(seq,error)
  best.idx <- which.min(error)
  best.dim <- M[best.idx]
  rf.fit <- randomForest(as.factor(class) ~ ., data=train.dat[,1:best.dim])
  #Ytrainhat = rf.fit$predicted
  preds <- predict(rf.fit,x.test)
  test.error <- mean(preds != y.test)
  return(test.error)
}

#Runs knearest for a given k
# Params: df: dataframe, effect: the classifier, kn = # neighbors
knearest <- function(df,effect,kn) {
  set.seed(23527)
  errors <- vector(mode='numeric', length=20)
  
  n <- dim(df)[1]
  
  for(j in 1:10){
    train.index <- sample(1:n, round(.9*n))
    x.train <- df[train.index,]
    y.train <- effect[train.index]
    train.dat <- data.frame(class = y.train, x.train)
    dim(train.dat)
    # testing set
    x.test <- df[-train.index,]
    y.test <- effect[-train.index]
    
    fit <- knn(train = x.train, test = x.test, cl = y.train, k = kn)
    errors[j] <- mean(fit != y.test)
  }
  
  return(mean(errors))
}

#First finds best k using CV then tests that k 
cv.knearest <- function(df,effect) {
  #Use cross-val to find k which minimizes CV error
  set.seed(23527)
  M <- seq(1,20,2)    # create the vector of tuning parameters
  n <- nrow(df)
  
  # create a vector to hold the 10 error rates
  cv.error <- vector(mode = "numeric", length = 10)
  
  # create a vector to hold the 10 best values of k
  best.k <- vector(mode = "numeric", length = 10)
  
  for(i in 1:10){
    
    train.ix <- sample(1:n, (9/10)*n)
    meth.train = df[train.ix,]
    effect.train = effect[train.ix]
    meth.test = df[-train.ix,]
    effect.test = effect[-train.ix]
    
    m.error <- vector(mode = "numeric")
    
    for(m in M){
      
      inside.error <- vector(mode = "numeric", length = 20)
      
      for(j in 1:20){
        n2 <- nrow(meth.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        meth.train2 = meth.train[train2.ix,]
        effect.train2 = effect.train[train2.ix]
        meth.validate = meth.train[-train2.ix,]
        effect.validate = effect.train[-train2.ix]
        
        fit <- knn(train = meth.train2, test = meth.validate, cl = effect.train2, k = m)
        inside.error[j] <- mean(fit != effect.validate)
      }
      
      m.error <- append(m.error, mean(inside.error))
    }
    
    # find the k with the smallest error
    idx <- which.min(m.error)
    best.k[i] <- M[idx]
    
    fit.train <- knn(train = meth.train, test = meth.validate, cl = effect.train, k = M[idx])
    cv.error[i] <- mean(fit.train != effect.test)
    
  }
  #run knn
  err_table <- cbind(best.k,cv.error)
  new_idx <- which.min(cv.error)
  top.k <- err_table[new_idx,'best.k']
  print(top.k)
  k.fit <- knn(train = meth.train, test = meth.test, cl = effect.train, k = top.k)
  test.error <- mean(fit.train != effect.test)
  return(test.error) #runs knn with the best k
}

 #runs svm using a linear separation kernel
svm_linear <- function(df,effect, costf = 1) {
  # split the centered data into training and testing
  set.seed(23527)
  errors <- vector(mode='numeric', length=10)
  
  n <- nrow(df)
  
  
  # training set
  for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- df[train.index,]
    y.train <- effect[train.index]
    train.dat <- data.frame(class = y.train, x.train)
    #dim(train.dat)
    # testing set
    x.test <- df[-train.index,]
    y.test <- effect[-train.index]
    test.dat <- data.frame(class = y.test, x.test)
    #dim(test.dat)
    svm_model <- svm(as.factor(y.train) ~ ., data=x.train, kernel="linear")
    #Ytrainhat = rf.fit$predicted
    preds <- predict(svm_model,x.test)
    #trainingerror$rf = 1 - mean(Ytrainhat==Ytrain)
    #testerror$rf = 1 - mean(Ytesthat==Ytest)
    
    # get the classification success rate
    errors[i] <- mean(preds != y.test)    
  }
  
  return(mean(errors))
}

#Does Cross validation to determine both the best number of features and the best value for cost function
cv.svm_linear <- function(df,effect) {
  set.seed(23527)
  dim_seq=seq(5,60,5)
  cost_seq=seq(0.25,2.5,0.25)
  
  # create the vector of tuning parameters
  n <- nrow(df)
  
  train.ix <- sample(1:n, (0.9)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
  cost.error <- matrix(0,nrow=12,ncol=10)
  
  for(i in 1:12){
    
    for(m in 1:10){
      
      inside.error <- vector(mode = "numeric", length = 5)
      #cv.features <- matrix(NA,nrow = 5,ncol = 100)
      
      for(j in 1:10){
        
        n2 <- nrow(Xtrain)
        train2.ix <- sample(1:n2, (0.9)*n2)
        Xtrain2 = Xtrain[train2.ix,]
        Ytrain2 = Ytrain[train2.ix]
        
        pvals <- data.frame(GENE = colnames(df))
        pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
        pvals_sel = pvals[order(pvals$p)[1:dim_seq[i]],]
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
  
  idx <- which(cost.error == min(cost.error), arr.ind = TRUE)
  best.dim <- dim_seq[idx[1]]
  best.cost <- cost_seq[idx[2]]
  return(svm_linear(df[,1:best.dim],effect,cost=best.cost))
}

#runs svm using a radial separation kernel
svm_rbf <- function(df,effect, cost=1) {
  # split the centered data into training and testing
  set.seed(23527)
  errors <- vector(mode='numeric', length=10)
  
  n <- nrow(pc.scores)
  
  
  # training set
  for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- df[train.index,]
    y.train <- effect[train.index]
    train.dat <- data.frame(class = y.train, x.train)
    #dim(train.dat)
    # testing set
    x.test <- df[-train.index,]
    y.test <- effect[-train.index]
    test.dat <- data.frame(class = y.test, x.test)
    #dim(test.dat)
    svm_model <- svm(as.factor(y.train) ~ ., data=x.train, kernel="radial")
    #Ytrainhat = rf.fit$predicted
    preds <- predict(svm_model,x.test)
    #trainingerror$rf = 1 - mean(Ytrainhat==Ytrain)
    #testerror$rf = 1 - mean(Ytesthat==Ytest)
    
    # get the classification success rate
    errors[i] <- mean(preds != y.test)    
  }
  
  return(mean(errors))
}

#Does Cross validation to determine both the best number of features and the best value for cost function
cv.svm_rbf <- function(df,effect) {
  set.seed(23527)
  dim_seq=Q
  cost_seq=seq(0.25,2.5,0.25)
  
  # create the vector of tuning parameters
  n <- nrow(df)
  
  train.ix <- sample(1:n, (0.9)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
  cost.error <- matrix(0,nrow=length(Q),ncol=10)
  
  for(i in 1:length(Q)){
    
    for(m in 1:10){
      
      inside.error <- vector(mode = "numeric", length = 5)
      #cv.features <- matrix(NA,nrow = 5,ncol = 100)
      
      for(j in 1:10){
        
        n2 <- nrow(Xtrain)
        train2.ix <- sample(1:n2, (0.9)*n2)
        Xtrain2 = Xtrain[train2.ix,]
        Ytrain2 = Ytrain[train2.ix]
        
        pvals <- data.frame(GENE = colnames(df))
        pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
        pvals_sel = pvals[order(pvals$p)[1:dim_seq[i]],]
        
        Xtrain2 = Xtrain2[,pvals_sel$GENE]
        Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
        Yvalidate = Ytrain[-train2.ix]
        
        svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="radial",cost=cost_seq[m])
        Ypred<-predict(svm_model,Xvalidate)
        inside.error[j] <- mean(Ypred != Yvalidate)
      }
      
      cost.error[i,m] = mean(inside.error) 
      
    }
  }
  
  idx <- which(cost.error == min(cost.error), arr.ind = TRUE)
  best.dim <- dim_seq[idx[1]]
  best.cost <- cost_seq[idx[2]]
  return(svm_rbf(df[,1:best.dim],effect,cost=best.cost))
}

# Jacob #
#### Naive Bayes, LASSO, Ridge, Elastic Net Classifiers:
# Returns best test error out of 10 trials
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

drugids <- c("1007", "1026", "1016", "1006", "1058", "1015", "1014", "1001", "1053",
             "1060", "1011", "1054", "1037", "1066", "1008")

# loading data
screening <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/screening.rds")
expr <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/expression.rds")

effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT")]
expr.1014 <- expr[as.character(effective.1014$CL),]

# creating expression data for each drug in drugids: 
# expr.{drugname} & effective.{drug}

# Lists to hold expr.{drugname} (sig. expression data for each drug) & 
# effective.{drug} (effectiveness for each CL w/ each drug)

list_effective <- list()
list_expr <- list()

for (i in 1:length(drugids)) {
  
  df_effective <- subset(screening, DRUG_ID_lib == drugids[i])[,c("CL","EFFECT")]
  df_expr <- expr[as.character(df_effective$CL),]
  
  list_effective[[i]] <- df_effective
  list_expr[[i]] <- df_expr
}

##### Performing tests for each drug for expression: to use, vary list_expr[[i]] for i in 1:15 for each drug...

# One data frame for t-test selection:
t_test_drug_expression <- data.frame(feature_selection_method = numeric(),
                                     number_of_features = numeric(),
                                     classifier = numeric(),
                                     parameters = numeric(),
                                     test_error = numeric())
pca_test_drug_expression <- data.frame(feature_selection_method = numeric(),
                                        number_of_features = numeric(),
                                        classifier = numeric(),
                                        parameters = numeric(),
                                        test_error = numeric())
## SELECTION OF DRUG: 

## t-test feature selection
pvals <- data.frame(GENE = colnames(list_expr[[1]]))
pvals$p = apply(list_expr[[1]],2,get.p, list_effective[[1]]$EFFECT)
pvals_sel = pvals[pvals$p<=.05,]
pvals_sel <- pvals_sel %>% arrange(p) 

# pca feature selection:
center.expr=scale(list_expr[[1]],scale=F)
pca <- svd(center.expr, nu = 20, nv = 0)  
m = sum(cumsum(pca$d^2)/sum(pca$d^2)<0.9)

pca_comp=prcomp(center.expr)
pc.scores <-as.data.frame(pca_comp$x)

# using genes we want (significant covariates)
list_expr[[1]] <- list_expr[[1]][,pvals_sel$GENE]

effect <- list_effective[[1]]$EFFECT
expr.drug <- list_expr[[1]]

Q <- seq(50, 1000, 50) # feature number
alpha <- seq(0, 1, 0.2) # alpha level

## T-test feature selection classifiers:
cv_naive_bayes(expr.drug, effect)
cv_lasso(expr.drug, effect)
cv_ridge(expr.drug, effect)
cv_elas(expr.drug, effect)

rand_forest(expr.drug, effect)
knearest(expr.drug, effect, 5)
svm_linear(expr.drug, effect)
svm_rbf(expr.drug, effect)

cv.knearest(expr.drug, effect)
cv.rand_forest(expr.drug, effect)
cv.svm_linear(expr.drug,effect)
cv.svm_rbf(expr.drug,effect)

## PCA feature selection classifiers: (use all features...)

p <- dim(pc.scores)[1]
Q <- seq(5, p, 5) # Feature number for pca
pc.scores

cv_naive_bayes(pc.scores, effect)
cv_lasso(as.data.frame(pc.scores), effect)
cv_ridge(pc.scores, effect)
cv_elas(pc.scores, effect)

rand_forest(pc.scores, effect)
knearest(pc.scores, effect, 5)
svm_linear(pc.scores, effect)
svm_rbf(pc.scores, effect)

cv.knearest(pc.scores, effect)
cv.rand_forest(pc.scores, effect, M = Q)
cv.svm_linear(pc.scores,effect)
cv.svm_rbf(pc.scores,effect)
















