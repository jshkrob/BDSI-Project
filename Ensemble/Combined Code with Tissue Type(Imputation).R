###############Combined using imputed copynumber data using missForest
library(matlib)
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
library(class)

library(tidyverse)
library(dplyr)
library(readxl)
library(varhandle)



pcaknn<-function(df,effect)
{
  
  library(caret)
  n=nrow(df)
  
  kinit=sample(seq(1,19,2),size=1)
  dim_seq=seq(5,n/5,5)
  kp_seq=seq(1,19,2)
  
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  Ytrain=effect[train.ix]
  length(Ytrain)
  dim(Xtrain)
  dim(Xtest)
  Ytest=effect[-train.ix]
  ##kp=c(kinit,0,0,0,0)
  kp=c(kinit,0,0,0,0)
  kp
  for(i in 1:4)##1:4
  {dim.error <- vector(mode = "numeric")
  kp.error <- vector(mode = "numeric")
  
  for(m in dim_seq)
  {
    inside.error <- vector(mode = "numeric", length = 5)
    for(j in 1:5)
    {
      n2 <- nrow(Xtrain)
      n2
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      Xvalidate=Xtrain[-train2.ix,]
      Yvalidate=Ytrain[-train2.ix]
      pccomp=prcomp(Xtrain2)
      pctrain=(pccomp$x)[,1:m]
      pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:m]
      train2 = data.frame(pctrain, Ytrain2)      # assemble dataframes
      test2 = data.frame(pctest, Yvalidate)
      
      model5a = knn3(Ytrain2 ~ ., data=train2,  k = kp[i]) 
      Ypred = ifelse(predict(model5a,data.frame(pctest))[,2]>=.5,1,0)
      
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim<- dim_seq[idx]
  
  for(q in kp_seq)
  {inside.error <- vector(mode = "numeric", length = 5)
  for(j in 1:5)
  {
    n2 <- nrow(Xtrain)
    train2.ix <- sample(1:n2, (4/5)*n2)
    Xtrain2 = Xtrain[train2.ix,]
    Ytrain2 = Ytrain[train2.ix]
    Xvalidate=Xtrain[-train2.ix,]
    Yvalidate=Ytrain[-train2.ix]
    pccomp=prcomp(Xtrain2)
    pctrain=(pccomp$x)[,1:best.dim]
    pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:best.dim]
    
    train2 = data.frame(pctrain, Ytrain2)      # assemble dataframes
    test2 = data.frame(pctest, Yvalidate)
    
    model5a = knn3(Ytrain2 ~ ., data=train2,  k = q) 
    Ypred = ifelse(predict(model5a,data.frame(pctest))[,2]>=.5,1,0)
    
    inside.error[j] <- mean(Ypred != Yvalidate)
    
  }
  kp.error=append(kp.error, mean(inside.error))
  }
  idx <- which.min(kp.error)
  best.kp<- kp_seq[idx]
  kp[i+1]=best.kp
  }
  
  best=c(best.kp,best.dim)
  best
  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[,1:best.dim]
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]
  train1 = data.frame(pctrain1, Ytrain)      # assemble dataframes
  test1 = data.frame(pctest1, Ytest)
  model5a = knn3(Ytrain ~ ., data=train1,  k = best.kp) 
  Ypred = ifelse(predict(model5a,data.frame(pctest1))[,2]>=.5,1,0)
  cv.error=mean(Ypred!=Ytest)
  cv.error
  return(c(best.kp,best.dim,cv.error))
  
}
pcanb<-function(df,effect)
{n=nrow(df)
dim_seq=seq(5,n/5,5)
library(naivebayes)

train.ix <- sample(1:n, 0.8*n)
Xtrain = df[train.ix,]
Ytrain =effect[train.ix]
Xtest = df[-train.ix,]
Ytest = effect[-train.ix]

dim.error <- vector(mode = "numeric")
for(m in dim_seq){
  
  inside.error <- vector(mode = "numeric", length = 5)
  
  for(j in 1:5){
    n2 <- nrow(Xtrain)
    train2.ix <- sample(1:n2, (4/5)*n2)
    Xtrain2 = Xtrain[train2.ix,]
    Ytrain2 = Ytrain[train2.ix]
    Xvalidate = Xtrain[-train2.ix,]
    Yvalidate = Ytrain[-train2.ix]
    
    
    pccomp=prcomp(Xtrain2)
    pctrain=pccomp$x
    pctest=as.matrix(Xvalidate)%*%pccomp$rotation
    model1 = naive_bayes(x=pctrain[,1:m],y =Ytrain2)
    Ypred = predict(model1,pctest[,1:m])
    inside.error[j] <- mean(Ypred != Yvalidate)
  }
  
  dim.error <- append(dim.error, mean(inside.error))
}
idx <- which.min(dim.error)
best.dim<- dim_seq[idx]
pccomp=prcomp(Xtrain)
pctrain=pccomp$x
pctest=as.matrix(Xtest)%*%pccomp$rotation


model2= naive_bayes(x=pctrain[,1:best.dim],y = Ytrain)
Ypred = predict(model2,pctest[,1:best.dim])
cv.error=mean(Ypred!= Ytest)

return(c(best.dim,cv.error))

}

pcarf<-function(df,effect)
{
  n=nrow(df)
  
  
  dim_seq=seq(5,n/5,5)
  library(randomForest)
  
  train.ix <- sample(1:n, 0.8*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest =effect[-train.ix]
  
  
  
  dim.error <- vector(mode = "numeric")
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      Xvalidate = Xtrain[-train2.ix,]
      Yvalidate = Ytrain[-train2.ix]
      
      
      pccomp=prcomp(Xtrain2)
      pctrain=pccomp$x
      pctest=as.matrix(Xvalidate)%*%pccomp$rotation
      train=data.frame(pctrain[,1:m],Ytrain2)
      test=data.frame(pctest[,1:m],Yvalidate)
      model6 <- randomForest(as.factor(Ytrain2) ~ ., data=train)
      Ypred = predict(model6,pctest[,1:m])
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim<- dim_seq[idx]
  pccomp=prcomp(Xtrain)
  pctrain=pccomp$x
  pctest=as.matrix(Xtest)%*%pccomp$rotation
  train=data.frame(pctrain[,1:best.dim],Ytrain)
  test=data.frame(pctest[,1:best.dim],Ytest)
  model6 <- randomForest(as.factor(Ytrain) ~ ., data=train)
  Ypred = predict(model6,pctest[,1:best.dim])
  
  cv.error=mean(Ypred!= Ytest)
  best.dim
  cv.error
  return(c(best.dim,cv.error))
  
}
pcalasso<-function(df,effect)
{
  library(glmnet)
  n=nrow(df)
  lambda_init=cv.glmnet(as.matrix(df),as.double(effect))$lambda.min
  lambda_init
  lambda=cv.glmnet(as.matrix(df),as.double(effect))$lambda
  lambda_seq=seq(min(lambda),max(lambda),length.out = 10)
  dim_seq=seq(5,n/5,5)
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  
  
  Ytrain=effect[train.ix]
  length(Ytrain)
  dim(Xtrain)
  dim(Xtest)
  Ytest=effect[-train.ix]
  lambda=c(lambda_init,0)
  lambda
  for(i in 1:2)
  {dim.error <- vector(mode = "numeric")
  lambda.error <- vector(mode = "numeric")
  
  for(m in dim_seq)
  {
    inside.error <- vector(mode = "numeric", length = 5)
    for(j in 1:5)
    {
      n2 <- nrow(Xtrain)
      n2
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      Xvalidate=Xtrain[-train2.ix,]
      Yvalidate=Ytrain[-train2.ix]
      
      
      pccomp=prcomp(Xtrain2)
      pctrain=(pccomp$x)[,1:m]
      pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:m]
      
      
      fit <- glmnet(x = as.matrix(pctrain), y = as.double(Ytrain2),family = "binomial",lambda = lambda[i])
      Ypred = predict(fit,as.matrix(pctest),type = "response")>.5
      length(Ypred)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim<- dim_seq[idx]
  
  for(q in lambda_seq)
  {inside.error <- vector(mode = "numeric", length = 5)
  for(j in 1:5)
  {
    n2 <- nrow(Xtrain)
    train2.ix <- sample(1:n2, (4/5)*n2)
    Xtrain2 = Xtrain[train2.ix,]
    Ytrain2 = Ytrain[train2.ix]
    Xvalidate=Xtrain[-train2.ix,]
    Yvalidate=Ytrain[-train2.ix]
    
    pccomp=prcomp(Xtrain2)
    pctrain=(pccomp$x)[,1:best.dim]
    pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:best.dim]
    
    fit <- glmnet(x = as.matrix(pctrain), y = as.double(Ytrain2),family = "binomial",lambda = q)
    Ypred = predict(fit,as.matrix(pctest),type = "response")>.5
    inside.error[j] <- mean(Ypred != Yvalidate)
    
  }
  lambda.error=append(lambda.error, mean(inside.error))
  }
  idx <- which.min(lambda.error)
  best.lambda<- lambda_seq[idx]
  lambda[2]=best.lambda
  }
  
  best=c(best.lambda,best.dim)
  best
  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[,1:best.dim]
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]
  
  fit <- glmnet(x = as.matrix(pctrain1), y = as.double(Ytrain),family = "binomial",lambda = best.lambda)
  Ypred = predict(fit,as.matrix(pctest1),type = "response")>.5
  cv.error=mean(Ypred!=Ytest)
  cv.error
  
  return(c(cv.error,best.dim,best.lambda))
  
}

pcaridge<-function(df,effect)
{ n=nrow(df)
library(glmnet)
dim_seq=seq(5,n/5,5)
lambda_init=cv.glmnet(as.matrix(df),as.double(effect),alpha=0)$lambda.min
lambda_init

lambda=cv.glmnet(as.matrix(df),as.double(effect),alpha=0)$lambda
lambda_seq=seq(min(lambda),max(lambda),length.out = 10)




train.ix <- sample(1:n, (4/5)*n)
Xtrain=df[train.ix,]
Xtest=df[-train.ix,]
Ytrain=effect[train.ix]
length(Ytrain)
dim(Xtrain)
dim(Xtest)
Ytest=effect[-train.ix]
lambda=c(lambda_init,0)
lambda
for(i in 1:2)
{dim.error <- vector(mode = "numeric")
lambda.error <- vector(mode = "numeric")

for(m in dim_seq)
{
  inside.error <- vector(mode = "numeric", length = 5)
  for(j in 1:5)
  {
    n2 <- nrow(Xtrain)
    n2
    train2.ix <- sample(1:n2, (4/5)*n2)
    Xtrain2 = Xtrain[train2.ix,]
    Ytrain2 = Ytrain[train2.ix]
    Xvalidate=Xtrain[-train2.ix,]
    Yvalidate=Ytrain[-train2.ix]
    
    
    pccomp=prcomp(Xtrain2)
    pctrain=(pccomp$x)[,1:m]
    pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:m]
    
    
    fit <- glmnet(x = as.matrix(pctrain), y = as.double(Ytrain2),family = "binomial",lambda = lambda[i],alpha=0)
    Ypred = predict(fit,as.matrix(pctest),type = "response")>.5
    length(Ypred)
    inside.error[j] <- mean(Ypred != Yvalidate)
    
  }
  dim.error <- append(dim.error, mean(inside.error))
}
idx <- which.min(dim.error)
best.dim<- dim_seq[idx]

for(q in lambda_seq)
{inside.error <- vector(mode = "numeric", length = 5)
for(j in 1:5)
{
  n2 <- nrow(Xtrain)
  train2.ix <- sample(1:n2, (4/5)*n2)
  Xtrain2 = Xtrain[train2.ix,]
  Ytrain2 = Ytrain[train2.ix]
  Xvalidate=Xtrain[-train2.ix,]
  Yvalidate=Ytrain[-train2.ix]
  pccomp=prcomp(Xtrain2)
  pctrain=(pccomp$x)[,1:best.dim]
  pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:best.dim]
  
  fit <- glmnet(x = as.matrix(pctrain), y = as.double(Ytrain2),family = "binomial",lambda = q,alpha=0)
  Ypred = predict(fit,as.matrix(pctest),type = "response")>.5
  inside.error[j] <- mean(Ypred != Yvalidate)
  
}
lambda.error=append(lambda.error, mean(inside.error))
}
idx <- which.min(lambda.error)
best.lambda<- lambda_seq[idx]
lambda[2]=best.lambda
}

best=c(best.lambda,best.dim)
best
pccomp=prcomp(Xtrain)
pctrain1=(pccomp$x)[,1:best.dim]
pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]

fit <- glmnet(x = as.matrix(pctrain1), y = as.double(Ytrain),family = "binomial",lambda = best.lambda,alpha=0)
Ypred = predict(fit,as.matrix(pctest1),type = "response")>.5
cv.error=mean(Ypred!=Ytest)
cv.error

return(c(cv.error,best.dim,best.lambda))

}

comblda<-function(df,effect)
{library(MASS)
  inside.error <- vector(mode = "numeric", length = 5)
  for( j in 1:5)
  {n=nrow(df)
  train.ix=sample(1:n,round(0.8*n))
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  model=lda(Xtrain,Ytrain)
  Ypred=predict(model,Xtest)$class
  inside.error[j]=mean(Ypred!=Ytest)
  }
  
  cv.errorlda=mean(inside.error)
  return(cv.errorlda)
}
combqda<-function(df,effect)
{library(MASS)
  inside.error <- vector(mode = "numeric", length = 5)
  for( j in 1:5)
  {n=nrow(data)
  train.ix=sample(1:n,round(0.8*n))
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  model=qda(Xtrain,Ytrain)
  Ypred=predict(model,Xtest)$class
  inside.error[j]=mean(Ypred!=Ytest)
  }
  
  cv.errorqda=mean(inside.error)
  return(cv.errorqda)
}
combnb<-function(df,effect)
{library(naivebayes)
  inside.error <- vector(mode = "numeric", length = 5)
  for( j in 1:5)
  {n=nrow(data)
  train.ix=sample(1:n,round(0.8*n))
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  model=naive_bayes(x=Xtrain,y=Ytrain)
  Ypred=predict(model,Xtest)
  inside.error[j]=mean(Ypred!=Ytest)              
  }
  
  cv.errornaivebayes=mean(inside.error)
  
  return(cv.errornaivebayes)
}
combrf<-function(df,effect)
{library(randomForest)
  inside.error <- vector(mode = "numeric", length = 5)
  for( j in 1:5)
  {n=nrow(data)
  train.ix=sample(1:n,round(0.8*n))
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
  train=data.frame(Xtrain,Ytrain)
  test=data.frame(Xtest,Ytest)
  model1=randomForest(as.factor(Ytrain)~.,data=train)
  Ypred=predict(model1,data.frame(Xtest))
  inside.error[j]=mean(Ypred!=Ytest)
  }
  
  cv.errorrf=mean(inside.error)
  return(cv.errorrf)
}
comblasso<-function(df,effect)
{library(glmnet)
  lambda = cv.glmnet(as.matrix(df),as.double(effect))$lambda
  lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
  n <- nrow(df)
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest =effect[-train.ix]
  
  lambda.error <- vector(mode = "numeric")
  
  for(m in lambda_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      Xvalidate = Xtrain[-train2.ix,]
      Yvalidate = Ytrain[-train2.ix]
      
      
      fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = m)
      Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error <- append(lambda.error, mean(inside.error))
  }
  
  # find the k with the smallest error
  idx <- which.min(lambda.error)
  best.lambda <- lambda_seq[idx]
  
  
  fit.train <- glmnet(x = as.matrix(Xtrain), y = as.double(Ytrain),family = "binomial",lambda = best.lambda )
  Ypred = predict(fit.train,as.matrix(Xtest),type = "response")>.5
  cv.error <- mean(Ypred!= Ytest)
  
  tablasso = cbind(best.lambda,cv.error) 
  return(tablasso)
}

combridge<-function(df,effect)
{library(glmnet)
  lambda = cv.glmnet(as.matrix(df),as.double(effect),alpha=0)$lambda
  lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
  n <- nrow(df)
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest =effect[-train.ix]
  
  lambda.error <- vector(mode = "numeric")
  
  for(m in lambda_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      Xvalidate = Xtrain[-train2.ix,]
      Yvalidate = Ytrain[-train2.ix]
      
      
      fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = m,alpha=0)
      Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error <- append(lambda.error, mean(inside.error))
  }
  
  # find the k with the smallest error
  idx <- which.min(lambda.error)
  best.lambda <- lambda_seq[idx]
  
  
  fit.train <- glmnet(x = as.matrix(Xtrain), y = as.double(Ytrain),family = "binomial",lambda = best.lambda,alpha=0 )
  Ypred = predict(fit.train,as.matrix(Xtest),type = "response")>.5
  cv.error <- mean(Ypred!= Ytest)
  
  tabridge = cbind(best.lambda,cv.error)
  return(tabridge)
}
combsvm<-function(df,effect)
{library(e1071)
  n=nrow(df)
  cost_seq=seq(0.25,2.5,0.25)
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest =effect[-train.ix]
  
  
  cost.error <- vector(mode = "numeric")
  for(m in cost_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      Xvalidate = Xtrain[-train2.ix,]
      Yvalidate = Ytrain[-train2.ix]
      
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=Xtrain2, kernel="linear",cost=m)
      Ypred = predict(svm_model,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    cost.error <- append(cost.error, mean(inside.error))
  }
  
  # find the k with the smallest error
  idx <- which.min(cost.error)
  best.cost <- cost_seq[idx]
  
  
  svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain, kernel="linear",cost=best.cost)
  Ypred = predict(svm_model,Xtest)
  cv.error <- mean(Ypred!= Ytest)
  
  tabsvm = cbind(best.cost,cv.error)
  return(tabsvm)
}

combknn<-function(df,effect)
{library(caret)
  n=nrow(df)
  n
  kp_seq=seq(1,19,2)
  kp.error=vector(mode = "numeric")
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest =effect[-train.ix]
  
  
  for(m in kp_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      
      Xvalidate = Xtrain[-train2.ix,]
      Yvalidate = Ytrain[-train2.ix]
      
      
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      model5a = knn3(Ytrain2 ~ ., data=train2,  k = m) 
      Ypred = ifelse(predict(model5a,data.frame(Xvalidate))[,2]>=.5,1,0)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    kp.error <- append(kp.error, mean(inside.error))
  }
  
  # find the k with the smallest error
  idx <- which.min(kp.error)
  best.kp <- kp_seq[idx]
  
  
  train = data.frame(Xtrain, Ytrain)      # assemble dataframes
  test= data.frame(Xtest, Ytest)
  model5a = knn3(Ytrain ~ ., data=train,  k = best.kp) 
  Ypred = ifelse(predict(model5a,data.frame(Xtest))[,2]>=.5,1,0)
  cv.error <- mean(Ypred!= Ytest)
  
  tabknn = cbind(best.kp,cv.error)
  return(tabknn)
}
######################

cell_line_type <- read_excel("Cell_Lines_Details.xlsx", sheet = "Cell line details" )
colnames(cell_line_type) <- c("CELL_LINE_NAME", "CL", "WES", "CNA", "GENE EXP", "METHY",
                              "Drug Response", "Site", "Tissue Site 2", "Cancer Type", "MSI", 
                              "Screen Medium", "Growth Properties")


cell_line_type <- cell_line_type[, c("CELL_LINE_NAME", "Site")]
screening <- readRDS("screening.rds")
expression <- readRDS("expression.rds")
meth <- readRDS("methylation_na.rds")
copynum <- readRDS("Mcopynum1014_0.1_KEEP_RMNA.impDATA.rds")
sum(is.na(meth))
screening <- screening %>% left_join(cell_line_type, by = "CELL_LINE_NAME")

# drug-specific data:
effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT", "Site","CELL_LINE_NAME")]
expr.1014 <- expression[as.character(effective.1014$CL), ]
meth.1014 <- meth[as.character(effective.1014$CELL_LINE_NAME), ]
Y <- effective.1014[, c("EFFECT", "Site")]
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
sum(is.na(meth.1014))
pvals <- data.frame(GENE = colnames(meth.1014))
pvals$p = apply(meth.1014,2,get.p,Y$EFFECT)
dim(meth.1014)
length(pvals$p)

pvals_sel = head(pvals[order(pvals$p),],10000)

meth.1014=meth.1014[,pvals_sel$GENE]
copynum.1014 <- data.frame(copynum[as.character(effective.1014$CL), ])
sum(is.na(copynum.1014))
## First, create a one-hot encoding for given drug: colnames = site.{tissue}
indicator_site <- data.frame(to.dummy(effective.1014$Site, prefix = "site"), 
                             row.names = effective.1014$CL)
## perform classifiers:


# choose either indicator_histology, indicator_site


# glm for "site" and "histography"


## Now subset (for the given drug) each data set (meth, expr, copynum) by the type of tissue, and create a "centered"
## data frame for each data set based on the tissue feature means:


expr.1014$Site <- effective.1014$Site # append Site column to group by tissues...
meth.1014$Site <- effective.1014$Site # append Site column to group by tissues...
copynum.1014$Site <- effective.1014$Site # append Site column to group by tissues...

scale_this <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  return(x - mu)
}

## Centering of data:
rm(meth)
# EXPR
expr.1014_c <- expr.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)
# METH
meth.1014_c <- meth.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)
##rownames(meth.1014_c) <- rownames(expr.1014)

# COPY
copynum.1014_c <- copynum.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)
##rownames(copynum.1014_c) <- rownames(expr.1014)

# accordingly, do the same with the effect (to remove infreqeunt tissue cases...)
sum(is.na(copynum.1014_c))

Y_c <- Y %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  ungroup() %>% 
  dplyr::select(-Site) %>% 
  pull(EFFECT)


indicator_site$Site=expr.1014$Site

indicator_site_c<-indicator_site %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  ungroup() %>% 
  dplyr::select(-Site)

## output: {dataset}.{drugid}_c and Y_c 
center.expr=scale(expr.1014_c,scale=F)
modelselecomb<-function(m,modelcen,modelind)
   {pccomp=prcomp(center.expr)
    pcexpr=pccomp$x[,1:m]
    a=pcexpr
    proj=a%*%Ginv(t(a)%*%a)%*%t(a)
    
    orthoproj=(diag(nrow(center.expr))-proj)
    orthometh=orthoproj%*%as.matrix(meth.1014_c)
    center.meth=scale(orthometh,scale=F)
    pccomp1=prcomp(center.meth)  
    pcmeth=pccomp1$x[,1:m]
    b=pcmeth
    c=cbind(a,b)
    
    proj1=c%*%Ginv(t(c)%*%c)%*%t(c)
    
    dim(proj1)
    orthoproj1=(diag(nrow(center.expr))-proj1)
    R=orthoproj1
    
    
    C1=R%*%as.matrix(copynum.1014_c)
    
    C1<-scale(C1,scale=F)
    
    sv=svd(C1,nu=m,nv=0)
    
    d=diag(sv$d[1:m])
    pccopy=sv$u%*%sqrt(d)
    data=cbind(pcexpr,pcmeth,pccopy)
    
    colnames(data)=as.character(c(1:(3*m)))
    colnames(data)
    if(modelcen=="lasso")
    {p=comblasso(data,Y_c)

    best.lambda=p[1]
    }else if(modelcen=="ridge")
    {p=combridge(data,Y_c)
    best.lambda=p[1]
    }else if (modelcen=="knn")
    { p=combknn(data,Y_c)
     best.kp=p[1]
    }
    
    if(modelind=="lasso")
    {
      lambda_ind=cv.glmnet(as.matrix(indicator_site_c),as.double(Y_c))$lambda.min
    }else(modelind=="ridge")
    {
      lambda_ind=cv.glmnet(as.matrix(indicator_site_c),as.double(Y_c))$lambda.min
    }
    
    centererror=rep(0,times=5)
    indicatorerror=rep(0,times=5)
    testerror=rep(0,times=5)
    for(i in 1:5)
    {n=nrow(center.expr)
    train.ix <- sample(1:n, (4/5)*n)
    Xtrain=data[train.ix,]
    Xtest=data[-train.ix,]
    Ytrain=Y_c[train.ix]
    Ytest=Y_c[-train.ix]
    centrain=data.frame(Xtrain,Ytrain)
    xtrain <- indicator_site_c[train.ix,]
    xtest <- indicator_site_c[-train.ix, ]
    indtrain=data.frame(xtrain,Ytrain)
    
    if(modelcen=="lasso")
    {fit <- glmnet(x = as.matrix(Xtrain), y = as.double(Ytrain),family = "binomial",lambda = best.lambda)
    Ypred = predict(fit,as.matrix(Xtest),type = "response")>.5
    centererror[i]=mean(Ypred!=Ytest)
    prob1train=predict(fit,as.matrix(Xtrain),type = "response")
    prob1test=predict(fit,as.matrix(Xtest),type = "response") 
    }else if(modelcen=="ridge")
    {
      fit <- glmnet(x = as.matrix(Xtrain), y = as.double(Ytrain),family = "binomial",lambda = best.lambda,alpha=0)
      Ypred = predict(fit,as.matrix(Xtest),type = "response")>.5
      centererror[i]=mean(Ypred!=Ytest)
      prob1train=predict(fit,as.matrix(Xtrain),type = "response")
      prob1test=predict(fit,as.matrix(Xtest),type = "response") 
    }else if(modelcen=="nb")
    {fitnb<-naiveBayes(Ytrain~.,data=centrain)
    Ypred=predict(fitnb,Xtest)
    centererror[i]=mean(Ypred!=Ytest)
    prob1train=predict(fitnb,Xtrain,type="raw")[,2]
    prob1test=predict(fitnb,Xtest,type="raw")[,2]
    }else if(modelcen=="rf")
    {
      model1=randomForest(as.factor(Ytrain)~.,data=centrain)
      Ypred=predict(model1,data.frame(Xtest))
      centererror[i]=mean(Ypred!=Ytest)
      prob1train=predict(model1,data.frame(Xtrain),type="prob")[,2]
      prob1test=predict(model1,data.frame(Xtest),type="prob")[,2]
    }else if(modelcen=="knn")
    {model5a = knn3(Ytrain ~ ., data=centrain,  k = best.kp) 
    Ypred=ifelse(predict(model5a,data.frame(Xtest))[,2]>=.5,1,0)
    centererror[i]=mean(Ypred!=Ytest)
    prob1train=predict(model5a,data.frame(Xtrain))[,2]
    prob1test=predict(model5a,data.frame(Xtest))[,2]
    }else if(modelcen=="lda")
    {model1 = lda(Xtrain,Ytrain)
    Ypred = predict(model1,Xtest)$class
    centererror[i]=mean(Ypred!=Ytest)
    prob1train=predict(model1,Xtrain)$posterior[,2]
    prob1test=predict(model1,Xtest)$posterior[,2]
    }
    
    
    if(modelind=="lasso")
    {fit <- glmnet(x = as.matrix(xtrain), y = as.double(Ytrain),family = "binomial",lambda = lambda_ind)
    Ypred = predict(fit,as.matrix(xtest),type = "response")>.5
    indicatorerror[i]=mean(Ypred!=Ytest)
    prob2train=predict(fit,as.matrix(xtrain),type = "response")
    prob2test=predict(fit,as.matrix(xtest),type = "response")
    }else if(modelind=="ridge")
    {fit <- glmnet(x = as.matrix(xtrain), y = as.double(Ytrain),family = "binomial",lambda = lambda_ind)
    Ypred = predict(fit,as.matrix(xtest),type = "response")>.5
    indicatorerror[i]=mean(Ypred!=Ytest)
    prob2train=predict(fit,as.matrix(xtrain),type = "response")
    prob2test=predict(fit,as.matrix(xtest),type = "response")
    }else if(modelind=="nb")
    {fitnb<-naiveBayes(Ytrain~.,data=indtrain)
    Ypred=predict(fitnb,xtest)
    indicatorerror[i]=mean(Ypred!=Ytest)
    prob2train=predict(fitnb,xtrain,type="raw")[,2]
    prob2test=predict(fitnb,xtest,type="raw")[,2]
    }else if(modelind=="rf")
    {model1=randomForest(as.factor(Ytrain)~.,data=indtrain)
    Ypred=predict(model1,data.frame(xtest))
    indicatorerror[i]=mean(Ypred!=Ytest)
    prob2train=predict(model1,data.frame(xtrain),type="prob")[,2]
    prob2test=predict(model1,data.frame(xtest),type="prob")[,2]
    }else if(modelind=="logistic")
    {model_glm <- glm(Ytrain ~ ., family = binomial(link = "logit"), data = indtrain,maxit=100)
    Ypred=ifelse(predict(model_glm, data.frame(xtest),type="response") >= 0.5, TRUE, FALSE)
    indicatorerror[i]=mean(Ypred!=Ytest)
    prob2train=predict(model_glm,data.frame(xtrain),type="response")
    prob2test=predict(model_glm,data.frame(xtest),type="response")
    }
    trainX=cbind(prob1train,prob2train)
    testX=cbind(prob1test,prob2test)
    colnames(testX)=c("1","2")
    colnames(trainX)=c("1","2")
    train=data.frame(trainX,Ytrain)
    test=data.frame(testX,Ytest)
    
    ##model=glmnet(x = as.matrix(trainX), y = as.double(Ytrain),family = "binomial")
    ##ytesthat   =predict(model,as.matrix(testX),type = "response")>.5
    #trainingerror = 1 - mean(ytrain==ytrainhat)
    model_glm <- glm(Ytrain ~ ., family = binomial(link = "logit"), data = train,maxit = 100)
    ytesthat <- ifelse(predict(model_glm, data.frame(testX),type="response") >= 0.5, TRUE, FALSE)
    testerror[i] =mean(Ytest!=ytesthat)
    coef1=model_glm$coefficients[2]
    coef2=model_glm$coefficients[3]
    }
    h=c(mean(coef1),mean(coef2),mean(testerror))
    return(h)
}
modelselecomb(10,"lda","logistic")
