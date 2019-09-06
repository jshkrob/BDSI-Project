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

# loading data
screening <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/screening.rds")
meth <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/methylation.rds")

# EDA: 
screening %>% 
  group_by(DRUG_ID_lib) %>% 
  summarise(size = n(), effective = sum(EFFECT), ineffective = sum(!EFFECT))

# packages for missing data: MICE, Amerlia, missForest, Hmisc, mi

### For now, examine drug 1014: 

effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT", "CELL_LINE_NAME")]
meth.1014 <- meth[as.character(effective.1014$CELL_LINE_NAME), ]

# remove columns with NA values
meth.1014_na <- meth.1014[ , !is.na(colSums(meth.1014))]
meth.1014 <- meth.1014_na

# feature selection:
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
pvals <- data.frame(METHSITE = colnames(meth.1014))
pvals$p = apply(meth.1014,2,get.p, effective.1014$EFFECT)
pvals_sel = pvals[pvals$p<=.05,]
pvals_best <- pvals %>% arrange(p) %>% slice(1:10000)
meth.1014 <- meth.1014[,pvals_best$METHSITE]

# pvals = pvals[complete.cases(pvals), ]

# Error data frame: repeat N (10) times test error and average for each method
ERROR <- data.frame(nb = numeric(10), lda = numeric(10), lasso = numeric(10), ridge = numeric(10), elas = numeric(10))

for(i in 1:10) {
  
  # choosing training and test
  n_obs <- nrow(meth.1014)
  train.ix = sample(1:n_obs, round(n_obs* 0.6))
  ytrain <- effective.1014$EFFECT[train.ix] 
  xtrain <- meth.1014[train.ix,]
  ytest <- effective.1014$EFFECT[-train.ix]
  xtest <- meth.1014[-train.ix, ]
  
  # naive bayes:
  model_nb = naive_bayes(x = xtrain, y = as.factor(ytrain))
  ytrainhat = predict(model_nb, xtrain)
  ytesthat = predict(model_nb, xtest)
  #train_error = 1 - mean(ytrain == ytrainhat)
  testerror = 1 - mean(ytest == ytesthat)
  ERROR[i,1] <- testerror  
  
  # lda:
  model_lda = lda(xtrain, ytrain)
  ytrainhat = predict(model_lda, xtrain)$class
  ytesthat = predict(model_lda, xtest)$class
  #trainingerror = 1 - mean(ytrain==ytrainhat)
  testerror = 1 - mean(ytest==ytesthat)
  ERROR[i,2] <- testerror  
  
  # # qda: using 50 features...
  # model_qda = qda(xtrain[,1:50], ytrain)
  # ytrainhat = predict(model_qda,xtrain[,1:50])$class
  # ytesthat = predict(model_qda,xtest[,1:50])$class
  # #trainingerror = 1 - mean(ytrain == ytrainhat)
  # testerror = 1 - mean(ytest == ytesthat)
  # ERROR[i,3] <- testerror  
  
  # # Logistic Regression: 25 features (10,000 is too much overfitting)
  # train <- data.frame(ytrain, xtrain[1:25])
  # model_glm <- glm(ytrain ~ ., family = binomial(link = "logit"), data = train, control = list(maxit = 100)) # increasing allowed iterations beyond default of 25...
  # ytrainhat <- ifelse(model_glm$fitted.values >= 0.5, 1, 0)
  # ytesthat <- ifelse(predict(model_glm, xtest) >= 0.5, TRUE, FALSE)
  # #trainingerror = 1 - mean(ytrain==ytrainhat)
  # testerror = 1 - mean(ytest==ytesthat)
  # ERROR[i,4] <- testerror  
  
  # Lasso
  model_lasso <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", alpha = 1)
  ytrainhat <- ifelse(predict(model_lasso, as.matrix(xtrain)) >= 0.5, TRUE, FALSE)
  ytesthat <- ifelse(predict(model_lasso, as.matrix(xtest)) >= 0.5, TRUE, FALSE)
  #trainingerror <- mean(ytrainhat != ytrain)
  testerror <- mean(ytesthat != ytest)
  ERROR[i,3] <- testerror  
  
  # Ridge
  model_ridge <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", alpha = 0)
  ytrainhat <- ifelse(predict(model_ridge, as.matrix(xtrain)) >= 0.5, TRUE, FALSE)
  ytesthat <- ifelse(predict(model_ridge, as.matrix(xtest)) >= 0.5, TRUE, FALSE)
  #trainingerror <- mean(ytrainhat != ytrain)
  testerror <- mean(ytesthat != ytest)
  ERROR[i,4] <- testerror  
  
  # Elastic Net:
  model_elas <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", alpha = 0.5)
  ytrainhat <- ifelse(predict(model_elas, as.matrix(xtrain)) >= 0.5, TRUE, FALSE)
  ytesthat <- ifelse(predict(model_elas, as.matrix(xtest)) >= 0.5, TRUE, FALSE)
  #trainingerror <- mean(ytrainhat != ytrain)
  testerror <- mean(ytesthat != ytest)
  ERROR[i,5] <- testerror
  
}

summary(ERROR)


#---------------------------------------------------------------Cross Validation----------------------------------------------------------

# feature number to use: 50-500 via p-values
Q <- seq(50,500,50)
n = nrow(meth.1014)

### Naive Bayes CV w/ feature number:
CV_NB <- data.frame(featurenum = numeric(20), error = numeric(20)) 

# 20 different observed best features w/ test error:
for (i in 1:20) {
  
  # training & test: 
  train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
  methylation.train = meth.1014[train.ix,]
  effective.train = effective.1014$EFFECT[train.ix]
  methylation.test = meth.1014[-train.ix,]
  effective.test = effective.1014$EFFECT[-train.ix]
  
  best_q <- 0 
  CV_error <- vector(mode = "numeric")
  for(q in Q) {
    
    cv_error <- vector(mode = "numeric", length = 20)
    for(j in 1:20) {
      
      # splitting data into k folds:
      n2 <- nrow(methylation.train)
      train2.ix <- sample(1:n2, (19/20)*n2)
      
      # 1 cv fold, 19 training folds:
      methylation.train2 = methylation.train[train2.ix,]
      effective.train2 = effective.train[train2.ix]
      methylation.validate = methylation.train[-train2.ix,]
      effective.validate = effective.train[-train2.ix]
      
      # Fit naive Bayes classifier on folds, compute CV error using validation set
      model_nb <- naive_bayes(x = methylation.train2[, 1:q], y = effective.train2)
      ytesthat <- predict(model_nb, methylation.validate[, 1:q])
      cv_error[j] <- 1 - mean(ytesthat == effective.validate)
      
    }
    # CV error for each feature number
    CV_error <- append(CV_error, mean(cv_error))
  }
  
  # Find q with the smallest error
  idx <- which.min(CV_error)
  CV_NB[i, ] <- c(Q[idx], CV_error[idx])
}

### Lasso CV w/ features and lamda
CV_LASSO <- data.frame(lambda = numeric(20), featurenum = numeric(20), error = numeric(20))
l_min <- 1

for(i in 1:20) {
  
  # for now, 6 iterations to find converging feature # and lamda:
  for(k in 1:6) {
    
    # FIRST, create TEST and TRAINING sets:
    n = nrow(meth.1014)
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    methylation.train = meth.1014[train.ix,]
    effective.train = effective.1014$EFFECT[train.ix]
    methylation.test = meth.1014[-train.ix,]
    effective.test = effective.1014$EFFECT[-train.ix]
    
    # SECOND, select best # of features (q) given lambda via CV error
    best_q <- 0
    best_error <- 1
    
    for(q in Q) { 
      # CV step:
      cv_error <- vector(mode = "numeric", length = 20)
      
      for(j in 1:20) {
        # splitting data into k folds:
        n2 <- nrow(methylation.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        methylation.train2 = methylation.train[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        
        methylation.validate = methylation.train[-train2.ix,]
        effective.validate = effective.train[-train2.ix]
        
        # train glm with regul. using q features w/ l_min choice & 19 training folds
        lasso.fit <- glmnet(x = as.matrix(methylation.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min)
        test_predic <- ifelse(predict(lasso.fit, as.matrix(methylation.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
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
    lasso.fit <- cv.glmnet(x = as.matrix(methylation.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 1)
    l_min <- lasso.fit$lambda.min
  }
  
  # use selected q and lamda on test set to estimtae error:
  lasso.fit <- glmnet(x = as.matrix(methylation.train[,1:best_q]), 
                      y = effective.train, 
                      family = "binomial", 
                      lambda = l_min)
  
  error <- mean((ifelse(predict(lasso.fit, as.matrix(methylation.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
  
  CV_LASSO[i,1] <- l_min
  CV_LASSO[i,2] <- best_q
  CV_LASSO[i,3] <- error
}

### Ridge CV w/ features and lamda
CV_RIDGE <- data.frame(lambda = numeric(20), featurenum = numeric(20), error = numeric(20))
l_min <- 1

for(i in 1:20) {
  
  # for now, 6 iterations to find converging feature # and lamda:
  for(k in 1:6) {
    
    # FIRST, create TEST and TRAINING sets:
    n = nrow(meth.1014)
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    methylation.train = meth.1014[train.ix,]
    effective.train = effective.1014$EFFECT[train.ix]
    methylation.test = meth.1014[-train.ix,]
    effective.test = effective.1014$EFFECT[-train.ix]
    
    # SECOND, select best # of features (q) given lambda via CV error
    best_q <- 0
    best_error <- 1
    
    for(q in Q) { 
      # CV step:
      cv_error <- vector(mode = "numeric", length = 20)
      
      for(j in 1:20) {
        # splitting data into k folds:
        n2 <- nrow(methylation.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        methylation.train2 = methylation.train[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        
        methylation.validate = methylation.train[-train2.ix,]
        effective.validate = effective.train[-train2.ix]
        
        # train ridge with regul. using q features w/ l_min choice & 19 training folds
        ridge.fit <- glmnet(x = as.matrix(methylation.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min, alpha = 0)
        test_predic <- ifelse(predict(ridge.fit, as.matrix(methylation.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
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
    ridge.fit <- cv.glmnet(x = as.matrix(methylation.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 0)
    l_min <- lasso.fit$lambda.min
  }
  
  # use selected q and lamda on test set to estimtae error:
  ridge.fit <- glmnet(x = as.matrix(methylation.train[,1:best_q]), y = effective.train, family = "binomial", lambda = l_min, alpha = 0)
  
  error <- mean((ifelse(predict(ridge.fit, as.matrix(methylation.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
  
  CV_LASSO[i,1] <- l_min
  CV_LASSO[i,2] <- best_q
  CV_LASSO[i,3] <- error
}





