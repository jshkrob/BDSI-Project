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

pcasvm<-function(df,effect)
{
  library(e1071)
  
  cost_init=1
  n=nrow(df)
  dim_seq=seq(5,n/5,5)
  cost_seq=seq(0.25,2.5,0.25)
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  Ytrain=effect[train.ix]
  length(Ytrain)
  dim(Xtrain)
  dim(Xtest)
  Ytest=effect[-train.ix]
  cost=c(1,0,0,0,0)
  ##cost=c(1,0)
  for(i in 1:4)##1:4
  {dim.error <- vector(mode = "numeric")
  cost.error <- vector(mode = "numeric")
  
  for(m in dim_seq)
  {
    inside.error <- vector(mode = "numeric", length = 5)
    for(j in 1:5)
    {
      n2 <- nrow(Xtrain)
      n2
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      dim(Xtrain2)
      
      Ytrain2 = Ytrain[train2.ix]
      length(Ytrain2)
      Xvalidate=Xtrain[-train2.ix,]
      dim(Xvalidate)
      Yvalidate=Ytrain[-train2.ix]
      length(Yvalidate)
      pccomp=prcomp(Xtrain2)
      pctrain=(pccomp$x)[,1:m]
      pctest=(as.matrix(Xvalidate)%*%pccomp$rotation)[,1:m]
      
      
      svm_model <- svm(as.factor(Ytrain2)~ ., data=pctrain, kernel="linear",cost=cost[i])##cost[i])
      Ypred = predict(svm_model,pctest)
      length(Ypred)
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim<- dim_seq[idx]
  
  for(q in cost_seq)
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
    
    svm_model <- svm(as.factor(Ytrain2)~ ., data=pctrain, kernel="linear",cost=q)
    Ypred = predict(svm_model,pctest)
    inside.error[j] <- mean(Ypred != Yvalidate)
    
  }
  cost.error=append(cost.error, mean(inside.error))
  }
  idx <- which.min(cost.error)
  best.cost<- cost_seq[idx]
  cost[i+1]=best.cost
  }
  
  best=c(best.cost,best.dim)
  best
  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[,1:best.dim]
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]
  
  svm_model <- svm(as.factor(Ytrain)~ ., data=pctrain1, kernel="linear",cost=best.cost)
  Ypred = predict(svm_model,pctest1)
  cv.error=mean(Ypred!=Ytest)
  cv.error
  
  
  tab_svm=c(cv.error,best.cost,best.dim)
  return(tab_svm)
}

pcalda<-function(df,effect)
{
  n=nrow(df)
  
  
  dim_seq=seq(5,n/5,5)
  library(MASS)
  
  train.ix <- sample(1:n, 0.8*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
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
      
      model1 = lda(pctrain[,1:m],Ytrain2)
      Ypred = predict(model1,pctest[,1:m])$class
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim<- dim_seq[idx]
  pccomp=prcomp(Xtrain)
  pctrain=pccomp$x
  pctest=as.matrix(Xtest)%*%pccomp$rotation
  model1 = lda(pctrain[,1:best.dim],Ytrain)
  Ypred = predict(model1,pctest[,1:best.dim])$class
  
  cv.error=mean(Ypred!= Ytest)
  
  cv.error
  
  return(c(best.dim,cv.error))
}
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

cell_line_type <- read_excel("Cell_Lines_Details.xlsx", sheet = "Cell line details" )
colnames(cell_line_type) <- c("CELL_LINE_NAME", "CL", "WES", "CNA", "GENE EXP", "METHY",
                              "Drug Response", "Site", "Tissue Site 2", "Cancer Type", "MSI", 
                              "Screen Medium", "Growth Properties")


cell_line_type <- cell_line_type[, c("CELL_LINE_NAME", "Site")]
screening <- readRDS("screening.rds")
expression <- readRDS("expression.rds")
meth <- readRDS("methylation_na.rds")
copynum <- readRDS("copynumber.rds")

screening <- screening %>% left_join(cell_line_type, by = "CELL_LINE_NAME")

# drug-specific data:
effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT", "Site","CELL_LINE_NAME")]
expr.1014 <- expression[as.character(effective.1014$CL), ]
meth.1014 <- meth[as.character(effective.1014$CELL_LINE_NAME), ]
copynum.1014 <- copynum[as.character(effective.1014$CL), ]

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
rownames(meth.1014_c) <- rownames(expr.1014)

# COPY
copynum.1014_c <- copynum.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)
rownames(copynum.1014_c) <- rownames(expr.1014)

# accordingly, do the same with the effect (to remove infreqeunt tissue cases...)

Y <- effective.1014[, c("EFFECT", "Site")]
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



######################################FUNCTION

modelselec<-function(modelcen,modelind)
{  if(modelcen=="lasso")
{p=pcalasso(center.expr,Y_c)
best.dim=p[2]
best.lambda=p[3]
}else if(modelcen=="ridge")
{p=pcaridge(center.expr,Y_c)
best.dim=p[2]
best.lambda=p[3]
}else if(modelcen=="nb")
{best.dim=pcanb(center.expr,Y_c)[1]
}else if (modelcen=="knn")
{ p=pcaknn(center.expr,Y_c)
best.dim=p[2]
best.kp=p[1]
}else if(modelcen=="rf")
{best.dim=pcarf(center.expr,Y_c)[1]
}else if(modelcen=="lda")
{best.dim=pcalda(center.expr,Y_c)[1]
}else if(modelcen=="svm")
{p=pcasvm(center.expr,Y_c)
best.dim=p[3]
best.cost=p[2]
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
  Xtrain=center.expr[train.ix,]
  Xtest=center.expr[-train.ix,]
  Ytrain=Y_c[train.ix]
  Ytest=Y_c[-train.ix]
  pccomp=prcomp(Xtrain)
  pctrain1=(pccomp$x)[,1:best.dim]
  centrain=data.frame(pctrain1,Ytrain)
  pctest1=(as.matrix(Xtest)%*%pccomp$rotation)[,1:best.dim]
  xtrain <- indicator_site_c[train.ix,]
  xtest <- indicator_site_c[-train.ix, ]
  indtrain=data.frame(xtrain,Ytrain)
  
  if(modelcen=="lasso")
  {fit <- glmnet(x = as.matrix(pctrain1), y = as.double(Ytrain),family = "binomial",lambda = best.lambda)
  Ypred = predict(fit,as.matrix(pctest1),type = "response")>.5
  centererror[i]=mean(Ypred!=Ytest)
  prob1train=predict(fit,as.matrix(pctrain1),type = "response")
  prob1test=predict(fit,as.matrix(pctest1),type = "response") 
  }else if(modelcen=="ridge")
  {
    fit <- glmnet(x = as.matrix(pctrain1), y = as.double(Ytrain),family = "binomial",lambda = best.lambda,alpha=0)
    Ypred = predict(fit,as.matrix(pctest1),type = "response")>.5
    centererror[i]=mean(Ypred!=Ytest)
    prob1train=predict(fit,as.matrix(pctrain1),type = "response")
    prob1test=predict(fit,as.matrix(pctest1),type = "response") 
  }else if(modelcen=="nb")
  {fitnb<-naiveBayes(Ytrain~.,data=centrain)
  Ypred=predict(fitnb,pctest1)
  centererror[i]=mean(Ypred!=Ytest)
  prob1train=predict(fitnb,pctrain1,type="raw")[,2]
  prob1test=predict(fitnb,pctest1,type="raw")[,2]
  }else if(modelcen=="rf")
  {
    model1=randomForest(as.factor(Ytrain)~.,data=centrain)
    Ypred=predict(model1,data.frame(pctest1))
    centererror[i]=mean(Ypred!=Ytest)
    prob1train=predict(model1,data.frame(pctrain1),type="prob")[,2]
    prob1test=predict(model1,data.frame(pctest1),type="prob")[,2]
  }else if(modelcen=="knn")
  {model5a = knn3(Ytrain ~ ., data=centrain,  k = best.kp) 
  Ypred=ifelse(predict(model5a,data.frame(pctest1))[,2]>=.5,1,0)
  centererror[i]=mean(Ypred!=Ytest)
  prob1train=predict(model5a,data.frame(pctrain1))[,2]
  prob1test=predict(model5a,data.frame(pctest1))[,2]
  }else if(modelcen=="lda")
  {model1 = lda(pctrain1,Ytrain)
   Ypred = predict(model1,pctest1)$class
   centererror[i]=mean(Ypred!=Ytest)
   prob1train=predict(model1,data.frame(pctrain1))$posterior[,2]
   prob1test=predict(model1,data.frame(pctest1))$posterior[,2]
  }else if(modelcen=="svm")
  {
    svm_model <- svm(as.factor(Ytrain)~ ., data=pctrain1, kernel="linear",cost=best.cost)
    Ypred = predict(svm_model,pctest1)
    centererror[i]=mean(Ypred!=Ytest)
    svm_model <- svm(as.factor(Ytrain)~ ., data=pctrain1, kernel="linear",cost=best.cost,probability=TRUE)
    prob1train=attr(predict(svm_model,pctrain1,probability=TRUE),"probabilites")
    prob1test=attr(predict(svm_model,pctest1,probability=TRUE),"probabilities")
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

modelselec("rf","lasso")


##########################
