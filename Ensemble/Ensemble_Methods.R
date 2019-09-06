# Packages:
library(pROC)
library(ROCR)
library(readr)
library(pls)
library(e1071)
library(HiDimDA)
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

# data:
screening <- readRDS("screening.rds")
expr <- readRDS("expression.rds")
meth <- readRDS("methylation.rds")
copy <- readRDS("copynumber.rds")

###-------------------------------------------------------------for specific drug, load expr, meth, copynum data ------------------------------------------------

## pvalue function
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

## drug id: 1008
effective.1008 <- subset(screening, DRUG_ID_lib == "1008")[,c("CL","EFFECT", "CELL_LINE_NAME")]
Y=effective.1008$EFFECT

### Expression

expr.1008 <- expr[as.character(effective.1008$CL), ]
pvals <- data.frame(GENE = colnames(expr.1008))
pvals$p = apply(expr.1008,2,get.p,Y)
pvals_sel = head(pvals[order(pvals$p),],10000)
expr = expr.1008[,pvals_sel$GENE]
center.expr <- scale(expr)

### Methylation

mds.1008 <- meth[as.character(effective.1008$CELL_LINE_NAME), ]
mds.1008 <- mds.1008[ , !is.na(colSums(meth))]
pvals <- data.frame(GENE = colnames(mds.1008))
pvals$p = apply(mds.1008,2,get.p,Y)
pvals_sel = head(pvals[order(pvals$p),],10000)
mds = mds.1008[,pvals_sel$GENE]
center.meth=scale(mds)

### Copy Number
copynum.1037 <- readRDS(...)

###-------------------------------------------------------------ensemble functions (given parameters, feature # from excel sheet) -------------------------------

################## for functions

## output: list containing (1) results, (2) predicted values, (3) predicted probabilities of effectiveness, (4) AUC, (5) F1 SCORE
## specific elements:

# results = "results"
# predicted = "predicted values"
# predfitted = "predicted probabilites of effectiveness"
# auc = "AUC"
# f1score = "F1 SCORE"
# specifiticities
# sensitivities

## to access elements, simply to the following: output_of_method$results, output_of_model$auc, etc....


df <- center.expr
n <- nrow(df)
INDX <- sample(1:n, (1/2)*n)

f1 <- function(predicted, expected) {
  predict <- ifelse(predicted >= 0.5, 1, 0)
  total_pos_predicted <- sum(predict)
  
  # Precision: fraction of positive predicted results that are actually positive
  precision <- sum(predict & expected) / total_pos_predicted
  
  # Recall: fraction of actual positive results that are predicted pos
  recall <- sum(predict & expected) / sum(expected)
  f1 <- (2*precision*recall) / (precision + recall)
}

pcasvm_ensemble<-function(df, effect, roc = F, stack = FALSE, train.ix, best.cost, best.dim) {
  library(e1071)
  n=nrow(df)
  dim_seq=seq(5,n/5,5)
  cost_seq=seq(0.25,2.5,0.25)
  
  Xtrain = df[train.ix,]
  Xtest = df[-train.ix,]
  Ytrain = effect[train.ix]
  Ytest = effect[-train.ix]
  
  pccomp = prcomp(Xtrain)
  pctrain1 = (pccomp$x)[,1:best.dim]
  pctest1 = (as.matrix(Xtest) %*% pccomp$rotation)[ ,1:best.dim]
  
  svm_model <- svm(as.factor(Ytrain) ~ ., data=pctrain1, kernel="linear",cost=best.cost)
  Ypred = predict(svm_model,pctest1)
  cv.error = mean(Ypred!=Ytest)
  
  if (roc == TRUE) {
    svm_model <- svm(as.factor(Ytrain) ~ ., data=pctrain1, kernel="linear", cost=best.cost, probability = TRUE )
    Ypredfitted <- predict(svm_model, pctest1, probability = TRUE)[,2]
    results <- c(best.cost,best.dim,cv.error)
    
    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    AUC <- roc()$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE, 
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  tab_svm=c(cv.error,best.cost,best.dim)
  return(tab_svm)
}
pcaknn_ensemble <-function(df, effect, roc = F, stack = FALSE, train.ix, best.kp, best.dim) {
  
  library(caret)
  n=nrow(df)
  
  kinit=sample(seq(1,19,2),size=1)
  dim_seq=seq(5,n/5,5)
  kp_seq=seq(1,19,2)

  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  Ytrain=effect[train.ix]
  
  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[,1:best.dim]
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]
  train1 = data.frame(pctrain1, Ytrain)      # assemble dataframes
  test1 = data.frame(pctest1, Ytest)
  model5a = knn3(Ytrain ~ ., data=train1,  k = best.kp) 
  Ypred = ifelse(predict(model5a,data.frame(pctest1))[,2]>=.5,1,0)
  cv.error=mean(Ypred!=Ytest)
  cv.error
  
  if(roc == TRUE) {
    Ypredfitted <- predict(model5a,data.frame(pctest1))[,2]
    results <- c(best.kp,best.dim,cv.error)

    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    plot(roc, col = "#fc8d59", add = TRUE)
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE,
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  return(c(best.kp,best.dim,cv.error))
  
}
pcanb_ensemble <-function(df, effect, roc = F, stack = FALSE, train.ix, best.dim){
  n=nrow(df)
  dim_seq=seq(5,n/5,5)
  library(naivebayes)
  
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]

  pccomp=prcomp(Xtrain)
  pctrain=pccomp$x
  pctest=as.matrix(Xtest)%*%pccomp$rotation
  
  
  model2= naive_bayes(x=pctrain[,1:best.dim],y = Ytrain)
  Ypred = predict(model2,pctest[,1:best.dim])
  cv.error=mean(Ypred!= Ytest)
  
  if(roc == TRUE) {
    
    Ypredfitted <- predict(model2, pctest[,1:best.dim], type = "prob")[,2]
    results <- c(best.dim,cv.error)

    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    plot(roc, col = "#e6ab02", add = TRUE)
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE, 
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  return(c(best.dim,cv.error))
  
}
pcalda_ensemble <-function(df, effect, roc = F, stack = FALSE,train.ix, best.dim) {
  n=nrow(df)

  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
  pccomp=prcomp(Xtrain)
  pctrain=pccomp$x
  pctest=as.matrix(Xtest)%*%pccomp$rotation
  model1 = lda(pctrain[,1:best.dim],Ytrain)
  Ypred = predict(model1,pctest[,1:best.dim])$class
  
  cv.error=mean(Ypred!= Ytest)
  
  if(roc == TRUE) {
    Ypredfitted <- predict(model1, pctest[,1:best.dim])$posterior[,2]
    results <- c(best.dim,cv.error)

    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    #pred <- prediction(as.numeric(Ypredfitted), Ytest)
    #perf <- performance(pred, "tpr", "fpr")
    
    if (stack == TRUE) {
      plot(roc, col = "grey82", ci = TRUE, add = TRUE)
      
    } else {
      plot(roc, col = "grey82")
      
    }
    
    # plot(x = 1- roc$specificities, (roc$sensitivities))
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE, 
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  return(c(best.dim,cv.error))
}
pcaqda_ensemble <-function(df, effect, roc = F, stack = FALSE, train.ix, best.dim) {
  n=nrow(df)
  
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
  pccomp=prcomp(Xtrain)
  pctrain=pccomp$x
  pctest=as.matrix(Xtest)%*%pccomp$rotation
  model1 = qda(pctrain[,1:best.dim],Ytrain)
  Ypred = predict(model1,pctest[,1:best.dim])$class
  
  cv.error=mean(Ypred!= Ytest)
  
  if(roc == TRUE) {
    Ypredfitted <- predict(model1, pctest[,1:best.dim])$posterior[,2]
    results <- c(best.dim,cv.error)
    
    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE, 
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  return(c(best.dim,cv.error))
}
pcarf_ensemble <-function(df, effect, roc = F, stack = FALSE, train.ix, best.dim) {
  n=nrow(df)
  library(randomForest)
  
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest =effect[-train.ix]

  pccomp=prcomp(Xtrain)
  pctrain=pccomp$x
  pctest=as.matrix(Xtest)%*%pccomp$rotation
  train=data.frame(pctrain[,1:best.dim],Ytrain)
  test=data.frame(pctest[,1:best.dim],Ytest)
  model6 <- randomForest(as.factor(Ytrain) ~ ., data=train)
  Ypred = predict(model6,pctest[,1:best.dim])
  
  cv.error=mean(Ypred!= Ytest)
  
  if (roc == TRUE) {
    Ypredfitted <- predict(model6, pctest[,1:best.dim], type = "prob")[,2]
    results <- c(best.dim,cv.error)

    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE,                  
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  return(c(best.dim,cv.error))
  
}
pcalasso_ensemble <-function(df, effect, roc = F, stack = FALSE, train.ix, best.lambda, best.dim) {
  library(glmnet)
  n=nrow(df)
 
  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  Ytrain=effect[train.ix]
  Ytest=effect[-train.ix]
 
  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[ ,1:best.dim]
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]
  
  fit <- glmnet(x = as.matrix(pctrain1), y = as.double(Ytrain), family = "binomial", lambda = best.lambda)
  Ypred = predict(fit, as.matrix(pctest1), type = "response") > .5
  cv.error=mean(Ypred != Ytest)
  
  if(roc == TRUE) {
    Ypredfitted <- predict(fit, as.matrix(pctest1), type = "response")
    results <- c(cv.error,best.dim,best.lambda)

    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
    
    # extractable information
    return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE, 
                specifities = roc$specificities, sensitivities = roc$sensitivities))
  }
  
  return(c(cv.error,best.dim,best.lambda))
  
}
pcaridge_ensemble <-function(df, effect, roc = F, stack = FALSE, train.ix, best.lambda, best.dim) { 
  n=nrow(df)
  library(glmnet)

  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  Ytrain=effect[train.ix]
  Ytest=effect[-train.ix]

  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[,1:best.dim]
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]

  fit <- glmnet(x = as.matrix(pctrain1), y = as.double(Ytrain),family = "binomial",lambda = best.lambda, alpha=0)
  Ypred = predict(fit,as.matrix(pctest1),type = "response") > .5
  cv.error=mean(Ypred!=Ytest)
  cv.error

  if(roc == TRUE) {
    Ypredfitted <- predict(fit, as.matrix(pctest1), type = "response")
    results <- c(cv.error,best.dim,best.lambda)

    # printing ROC curve for specific model:
    roc <- roc(Ytest, as.numeric(Ypredfitted))
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
    
    AUC <- roc$auc
    F1_SCORE <- f1(Ypredfitted, Ytest)
  
  # extractable information
  return(list(results = cv.error, predicted = Ypred, predfitted = Ypredfitted, auc = AUC, f1score = F1_SCORE,
              specifities = roc$specificities, sensitivities = roc$sensitivities))
}

return(c(cv.error,best.dim,best.lambda))

}

# ensemble vote method: Chris
ensemble_vote <- function(effect, roc = F, stack = FALSE, train.ix, meth.fit, expr.fit, cn.fit) {
  Ytest <- effect[-train.ix]
  predic_ensemble <- data.frame(ifelse(meth_fit$predicted == TRUE, 1, 0),
                                ifelse(expr_fit$predicted == TRUE, 1, 0),
                                ifelse(copy_fit$predicted == TRUE, 1, 0))
  colnames(predic_ensemble) <- c("METH", "EXPR", "COPY")
  
  Ypred_ensemble <- predic_ensemble %>% 
    mutate_all(funs(sign_chg(.))) %>% 
    mutate(vote = is.positive(METH + EXPR + COPY)) %>% 
    pull(vote)
  
  ensemble_error <- mean(Ypred_ensemble != Ytest)
  
  roc <- roc(Ytest, as.numeric(Ypred_ensemble))
  
  if (roc) {
    
    if (stack == TRUE) {
      plot(roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(roc, col = "#d53e4f")
      
    }
  }
  
  AUC <- roc$auc
  F1_SCORE <- f1(Ypred_ensemble, Ytest)
  cv.error=mean(Ypred_ensemble != Ytest)
  sensitivities <- roc$sensitivies
  specificities <- roc$specificities
  specificities <-  1 - specificities
  return(list(results = cv.error, predicted = Ypred_ensemble, auc = AUC, f1score = F1_SCORE, 
              sensitivies = sensitivities, specificities=specificities))
}

#########################################

###---------------------------------------------------------tests------------------------------------------------------------------------------------------------
 
### Beforehand, load the methods that you want ahead of time

n <- nrow(center.expr)

spec_1 <- vector(mode = "numeric")
sens_1 <- vector(mode = "numeric")
for (i in 1:10) {
  n <- nrow(df)
  INDX <- sample(1:n, (4/5)*n)
  if (i == 1) {
    v <- pcalda_ensemble(center.expr, Y, train.ix = INDX, roc = TRUE, best.dim = 10) # load method 1 here
    spec_1 <- v$specifities
    sens_1 <- v$sensitivities
  } else {
    v <- pcalda_ensemble(center.expr, Y, train.ix = INDX, roc = TRUE, best.dim = 10) # load method 1 here
    spec_1 <- spec_1 + v$specifities
    sens_1 <- sens_1 + v$sensitivities
  }
}
avg_spec_1 <- spec_1/10
avg_sens_1 <- sens_1 / 10

#plot(1 - avg_spec_1, avg_sens_1)
df_1 <- data.frame(FPR = (1-avg_spec_1), TPR = avg_sens_1)


spec_2 <- vector(mode = "numeric")
sens_2 <- vector(mode = "numeric")
for (i in 1:10) {
  INDX <- sample(1:n, (4/5)*n)
  if (i == 1) {
    v <- pcalda_ensemble(center.expr, Y, train.ix = INDX, roc = TRUE, best.dim = 20) # load method here
    
    spec_2 <- v$specifities
    sens_2 <- v$sensitivities
  } else {
    v <- pcalda_ensemble(center.expr, Y, train.ix = INDX, roc = TRUE, best.dim = 20) # load method here
    
    spec_2 <- spec_2 + v$specifities
    sens_2 <- sens_2 + v$sensitivities
  }
}
avg_spec_2 <- spec_2/10
avg_sens_2 <- sens_2 / 10

#plot(1 - avg_spec_2, avg_sens_2)
df_2 <- data.frame(FPR = (1-avg_spec_2), TPR = avg_sens_2)

spec_3 <- vector(mode = "numeric")
sens_3 <- vector(mode = "numeric")
for (i in 1:10) {
  INDX <- sample(1:n, (4/5)*n)
  if (i == 1) {
    v <- pcalda_ensemble(center.expr, Y, train.ix = INDX, roc = TRUE, best.dim = 30) # load method here
    
    spec_3 <- v$specifities
    sens_3 <- v$sensitivities
  } else {
    v <- pcalda_ensemble(center.expr, Y, train.ix = INDX, roc = TRUE, best.dim = 30) # load method here
    
    spec_3 <- spec_3 + v$specifities
    sens_3 <- sens_3 + v$sensitivities
  }
}
avg_spec_3 <- spec_3/10
avg_sens_3 <- sens_3 / 10

#plot(1 - avg_spec_3, avg_sens_3)
df_3 <- data.frame(FPR = (1-avg_spec_3), TPR = avg_sens_3)

ggplot() +
  geom_line(data = df_1, aes(x = FPR, y = TPR), color = "blue", size = 1.0, alpha = 0.7) + 
  geom_line(data = df_2, aes(x = FPR, y = TPR), color = "green", size = 1.0, alpha = 0.7) + 
  geom_line(data = df_3, aes(x = FPR, y = TPR), color = "red", size = 1.0, alpha = 0.7) + 
  labs(title = "ROC Curve")

###------------------------------------------------------ensemble method: voted majority of result of each classifer ---------------------------------------------

meth_fit <- pcalasso_ensemble(center.meth, Y, roc = TRUE, stack = FALSE, train.ix = INDX, best.lambda = 0.03126065, best.dim = 25)
expr_fit <- pcalasso_ensemble(center.expr, Y, roc = TRUE, stack = TRUE, train.ix = INDX, best.lambda = 0.08327653, best.dim = 30)
copy_fit <- pcarf_ensemble(center.copy, Y, roc = TRUE, stack = TRUE, train.ix = INDX, best.dim = 5)

# make a legend:


# data frame of predictions
Ytest <- Y[-INDX]

predic_ensemble <- data.frame(ifelse(meth_fit$predicted == TRUE, 1, 0),
                              ifelse(expr_fit$predicted == TRUE, 1, 0),
                              ifelse(copy_fit$predicted == TRUE, 1, 0))
colnames(predic_ensemble) <- c("METH", "EXPR", "COPY")

# changes predictions of 0 to predictions of -1
sign_chg <- function(x) {
  ifelse(x == 0, -1, 1)
}
is.positive <- function(x) {
  x > 0
}

Ypred_ensemble <- predic_ensemble %>% 
  mutate_all(funs(sign_chg(.))) %>% 
  mutate(vote = is.positive(METH + EXPR + COPY)) %>% 
  pull(vote)

ensemble_error <- mean(Ypred_ensemble != Ytest)

# ensemble error via voting: 
ensemble_error


###------------------------------------------------------ensemble method: weighted average of result of each classifer (using 1-error)---------------------------------------------




###------------------------------------------------------ensemble method: baseline indicator matrix + centered data matrix ensemble---------------------------------------------

