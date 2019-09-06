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

#----------------------------------------------------PC Scores---------------------------------------
# first, we need to center the expression data by column means
meth <- meth.1014

center.meth = scale(meth, scale = F) #scale = F centers the data
# now use svd to perform PCA
pca <- svd(center.meth, nu = 56, nv = 0)
hist(log(pca$d))## plots log of distribution of sigular values (variance of the singular components)

m = sum(cumsum(pca$d^2)/sum(pca$d^2)<.8) #Looks at the cumulative sum vector of the variances over the sum of variance and finds the number of components to explain 80% of the variation
pc.scores <- as.data.frame(diag(pca$d)%*%pca$u) #computes pc score
pc.scores <- cbind(effective.1014, pc.scores) #matches to labels of data
m
dim(pc.scores)

head(pc.scores)

pca_comp = prcomp(center.meth)
pc.scores <- pca_comp$x[,1:60]
pc.scores <- data.frame(pc.scores)
#----------------------------------------------------Cross Validation----------------------------------------------------------

# feature number to use: 50-500 via p-values
Q <- seq(50,500,50)
alpha <- seq(0,1,0.2)
n = nrow(meth.1014)
  
  ### Naive Bayes CV w/ feature number:
CV_NB <- data.frame(featurenum = numeric(10), error = numeric(10)) 
for (i in 1:10) {
    
    # training & test: 
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    methylation.train = meth.1014[train.ix,]
    effective.train = effective.1014$EFFECT[train.ix]
    methylation.test = meth.1014[-train.ix,]
    effective.test = effective.1014$EFFECT[-train.ix]
    
    best_q <- 0 
    CV_error <- vector(mode = "numeric")
    for(q in Q) {
      
      cv_error <- vector(mode = "numeric", length = 10)
      for(j in 1:10) {
        
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
    best_q <- Q[idx]
    
    # Fit training data on selected feature number
    model_nb <- naive_bayes(x = methylation.train[, 1:best_q], y = effective.train)
    ytesthat <- predict(model_nb, methylation.test[, 1:q])
    test_error <- 1 - mean(ytesthat == effective.test)
    
    CV_NB[i, ] <- c(best_q, test_error)
  }
  
  ### Lasso CV 
CV_LASSO <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
l_min <- 1
for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:3) {
      
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
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
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
  
  ### Ridge CV 
CV_RIDGE <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
l_min <- 1
for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:3) {
      
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
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
          # splitting data into k folds:
          n2 <- nrow(methylation.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          methylation.train2 = methylation.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          
          methylation.validate = methylation.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # train glm with regul. using q features w/ l_min choice & 19 training folds
          lasso.fit <- glmnet(x = as.matrix(methylation.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min, alpha = 0)
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
      lasso.fit <- cv.glmnet(x = as.matrix(methylation.train[,1:best_q]), y = effective.train, family = "binomial", nfolds = 10, alpha = 0)
      l_min <- lasso.fit$lambda.min
    }
    
    # use selected q and lamda on test set to estimtae error:
    lasso.fit <- glmnet(x = as.matrix(methylation.train[,1:best_q]), 
                        y = effective.train, 
                        family = "binomial", 
                        lambda = l_min)
    
    error <- mean((ifelse(predict(lasso.fit, as.matrix(methylation.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_RIDGE[i,1] <- l_min
    CV_RIDGE[i,2] <- best_q
    CV_RIDGE[i,3] <- error
  }
  
  ### Elastic CV 
CV_ELAS <- data.frame(alpha = numeric(10), lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
l_min <- 1
  
for(i in 1:10) {
    
      # FIRST, create TEST and TRAINING sets:
      n = nrow(meth.1014)
      train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
      methylation.train = meth.1014[train.ix,]
      effective.train = effective.1014$EFFECT[train.ix]
      methylation.test = meth.1014[-train.ix,]
      effective.test = effective.1014$EFFECT[-train.ix]
      
      MCV <- vector(mode = "numeric")
      LAMBDA <- vector(mode = "numeric")
      
      # SECOND, select alpha and lambda based on CV:
      for(a in alpha) {
        # 10 fold cross validation for alpha = 0.0, 0.1, ..., 0.9, 1
        model_elas <- cv.glmnet(x = as.matrix(methylation.train), 
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
          n2 <- nrow(methylation.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          methylation.train2 = methylation.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          methylation.validate = methylation.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # Train using q features w/ l_min choice & 19 training folds
          elas.fit <- glmnet(x = as.matrix(methylation.train2[,1:q]), y = effective.train2, 
                              family = "binomial", lambda = lambda_best, alpha = alpha_best)
          
          test_predic <- ifelse(predict(elas.fit, as.matrix(methylation.validate[,1:q]),type = "response") >= 0.5, TRUE, FALSE)
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
      elas.fit <- glmnet(x = as.matrix(methylation.train[,1:best_q]), y = effective.train, 
                          family = "binomial", lambda = lambda_best, alpha = alpha_best)
      
      error <- mean((ifelse(predict(elas.fit, as.matrix(methylation.test[1:best_q])) >= 0.5, TRUE, FALSE)) != effective.test)
      
      CV_ELAS[i,1] <- alpha_best
      CV_ELAS[i,2] <- lambda_best
      CV_ELAS[i,3] <- best_q
      CV_ELAS[i,4] <- error
    }
    
  
 #----------------------------------------------------------------------Non-PCA Functions-------------------------------------------

## functions:
Q <- seq(50, 1000, 50)
cv_naive_bayes <- function(df, effect) {
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
    
    return(BEST_ERROR)
  }
cv_lasso <- function(df, effect) {
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
    print(CV_LASSO)
    return(BEST_ERROR)
    
  }
cv_ridge <- function(df, effect) {
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
    return(BEST_ERROR)
    
  }
cv_elas <- function(df, effect) {
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
    return(BEST_ERROR)
    
    
  }
  
cv_naive_bayes(meth.1014, effective.1014$EFFECT)
cv_lasso(meth.1014, effective.1014$EFFECT)
cv_ridge(meth.1014, effective.1014$EFFECT)
cv_elas(meth.1014, effective.1014$EFFECT)

pc.scores <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Methylation/meth_pcscores.rds")

Q <- seq(5,60,5)
cv_naive_bayes(pc.scores, effective.1014$EFFECT)
cv_lasso(pc.scores, effective.1014$EFFECT)
cv_ridge(pc.scores, effective.1014$EFFECT)
cv_elas(pc.scores, effective.1014$EFFECT)



























#--------------------------------------------------------------Using PCA/PCR-------------------------------------------------------

# These functions do not find an optimal number of featuers; instead, they use 
# the PCs for regression:

# Returns: CV error (TRAINING data) and test error (TEST data)

# Selects NUMB of features using CV, outputs best test error
cv.naive_bayes <- function(df, effect) {
  n = nrow(df)
  
  # training & test: 
  train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
  df.train = df[train.ix,]
  effective.train = effect[train.ix]
  df.test = df[-train.ix,]
  effective.test = effect[-train.ix]
  
  # cross validation error
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
    model_nb <- naive_bayes(x = df.train2, y = as.factor(effective.train2))
    ytesthat <- predict(model_nb, df.validate)
    cv_error[j] <- 1 - mean(ytesthat == effective.validate)
    
  }
  # mean CV error for each feature number
  CV_error <- append(CV_error, mean(cv_error))
  
  # fit training data and calculate TEST error
  model_nb <- naive_bayes(x = df.train, y = effective.train)
  ytesthat <- predict(model_nb, df.test)
  test_error <- 1 - mean(ytesthat == effective.test)
  
  # Return CV and test error for df
  return(c(CV_error,test_error))
}

cv.lasso <- function(df, effect) {
  CV_LASSO <- data.frame(lambda = numeric(10), error = numeric(10))
  l_min <- 1
  
  # 10 estimates of error 
  for(i in 1:10) {
    # FIRST, create TEST and TRAINING sets:
    n = nrow(df)
    train.ix = sample(1:n, round(n*0.60)) ## consider a approx. 60-40 split
    df.train = df[train.ix,]
    effective.train = effect[train.ix]
    df.test = df[-train.ix,]
    effective.test = effect[-train.ix]
   
    # Use feature number q to find best lamda for that q: 
    lasso.fit <- cv.glmnet(x = as.matrix(df.train), y = effective.train, family = "binomial", nfolds = 10, alpha = 1)
    l_min <- lasso.fit$lambda.min
  
    
    # Use selected q and lamda on test set to estimtae error:
    lasso.fit <- glmnet(x = as.matrix(df.train), y = effective.train, family = "binomial", lambda = l_min)
    
    error <- mean((ifelse(predict(lasso.fit, as.matrix(df.test)) >= 0.5, TRUE, FALSE)) != effective.test)
    
    CV_LASSO[i,1] <- l_min
    CV_LASSO[i,2] <- error
  }
  
  # return smallest error found among 10 examples:
  BEST_ERROR <- min(CV_LASSO$error)
  print(CV_LASSO)
  return(BEST_ERROR)
}

cv.ridge <- function(df, effect) {
  ### Ridge CV 
  CV_RIDGE <- data.frame(lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  for(i in 1:10) {
    
    # for now, 6 iterations to find converging feature # and lamda:
    for(k in 1:3) {
      
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
        cv_error <- vector(mode = "numeric", length = 10)
        
        for(j in 1:10) {
          # splitting data into k folds:
          n2 <- nrow(methylation.train)
          train2.ix <- sample(1:n2, (19/20)*n2)
          
          # 1 cv fold, 19 training folds:
          methylation.train2 = methylation.train[train2.ix,]
          effective.train2 = effective.train[train2.ix]
          
          methylation.validate = methylation.train[-train2.ix,]
          effective.validate = effective.train[-train2.ix]
          
          # train glm with regul. using q features w/ l_min choice & 19 training folds
          lasso.fit <- glmnet(x = as.matrix(methylation.train2[,1:q]), y = effective.train2, family = "binomial", lambda = l_min, alpha = 1)
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
    
    CV_RIDGE[i,1] <- l_min
    CV_RIDGE[i,2] <- best_q
    CV_RIDGE[i,3] <- error
  }
  # return smallest error found among 10 examples:
  BEST_ERROR <- min(CV_RIDGE$error)
  return(BEST_ERROR)
  
}

cv.elastic <- function(df, effect) {
  CV_ELAS <- data.frame(alpha = numeric(10), lambda = numeric(10), featurenum = numeric(10), error = numeric(10))
  l_min <- 1
  
  for(i in 1:10) {
    
    # FIRST, create TEST and TRAINING sets:
    n = nrow(meth.1014)
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
        n2 <- nrow(methylation.train)
        train2.ix <- sample(1:n2, (19/20)*n2)
        
        # 1 cv fold, 19 training folds:
        df.train2 = df[train2.ix,]
        effective.train2 = effective.train[train2.ix]
        df.validate = methylation.train[-train2.ix,]
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
  return(BEST_ERROR)
  
}

#############################################################################################

# PC Scores:

#####











