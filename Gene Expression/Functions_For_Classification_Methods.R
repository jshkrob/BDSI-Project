screening<-readRDS("screening.rds")
expr=readRDS("expression.rds")
effective.1014<- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT")]
expr.1014 <- expr[as.character(effective.1014$CL),]
Y=effective.1014$EFFECT
Y
center.expr=scale(expr.1014,scale=F)



###########Function for pca and svm using iterative method(Here the df should be centered)
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
cost=c(1,0)
cost
for(i in 1:2)
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
cost[2]=best.cost
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
return(c(cv.error,best.cost,best.dim))
}

pcasvm(center.expr,Y)

############################################################pcaknn(iterative)

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
  kp=c(kinit,0)
  kp
  for(i in 1:2)
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
  kp[2]=best.kp
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
pcaknn(center.expr,Y)
#########################################################pcanaivebayes


pcanb<-function(df,effect){
n=nrow(df)
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
pcanb(center.expr,Y)
###########################pcalda

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

pcalda(center.expr,Y)
###############################pcaqda

pcaqda<-function(df,effect)
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
      
      model1 = qda(pctrain[,1:m],Ytrain2)
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
  model1 = qda(pctrain[,1:best.dim],Ytrain)
  Ypred = predict(model1,pctest[,1:best.dim])$class
  
  cv.error=mean(Ypred!= Ytest)
  
  cv.error
  
  return(c(best.dim,cv.error))
}
pcaqda(center.expr,Y)
####################################pcarf

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
pcarf(center.expr,Y)


#########################logistic_glmnet_lasso_pca

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
  Ytest=Y[-train.ix]
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

pcalasso(center.expr,Y)
#########################logistic_glmnet_ridge_pca

pcaridge<-function(df,effect)
{
  library(glmnet)
  dim_seq=seq(5,n/5,5)
  lambda_init=cv.glmnet(as.matrix(df),as.double(effect),alpha=0)$lambda.min
  lambda_init
  
  lambda=cv.glmnet(as.matrix(df),as.double(effect),alpha=0)$lambda
  lambda_seq=seq(min(lambda),max(lambda),length.out = 10)
  
  
  n=nrow(df)
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain=df[train.ix,]
  Xtest=df[-train.ix,]
  Ytrain=effect[train.ix]
  length(Ytrain)
  dim(Xtrain)
  dim(Xtest)
  Ytest=Y[-train.ix]
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

##############################################################Pvalue(data is not centered.We take df as the original data.)

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


################################t-test naivebayes

nbtst<-function(df,effect)
{
  library(naivebayes)
  n=nrow(df)
  cv.error <- vector(mode = "numeric", length = 5)
  best.dim <- vector(mode = "numeric", length = 5)
  dim_seq=seq(50,500,50)
  dim_seq
  for(i in 1:5){
    
    train.ix <- sample(1:n, (4/5)*n)
    Xtrain = df[train.ix,]
    Ytrain = effect[train.ix]
    Xtest = df[-train.ix,]
    Ytest = effect[-train.ix]
    
    dim.error <- vector(mode = "numeric")
    
    for(m in dim_seq){
      
      inside.error <- vector(mode = "numeric", length = 5)
      #cv.features <- matrix(NA,nrow = 5,ncol = 100)
      
      for(j in 1:5){
        n2 <- nrow(Xtrain)
        train2.ix <- sample(1:n2, (4/5)*n2)
        Xtrain2 = Xtrain[train2.ix,]
        Ytrain2 = Ytrain[train2.ix]
        
        pvals <- data.frame(GENE = colnames(df))
        pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
        pvals_sel = head(pvals[order(pvals$p),],m)
        #cv.features[j,] = pvals_sel$GENE
        Xtrain2 = Xtrain2[,pvals_sel$GENE]
        Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
        Yvalidate = Ytrain[-train2.ix]
        
        model1 = naive_bayes(x=Xtrain2,y =Ytrain2)
        Ypred = predict(model1,Xvalidate)
        inside.error[j] <- mean(Ypred != Yvalidate)
      }
      
      dim.error <- append(dim.error, mean(inside.error))
    }
    idx <- which.min(dim.error)
    best.dim[i] <- dim_seq[idx]
    pvals <- data.frame(GENE = colnames(df))
    pvals$p = apply(Xtrain,2,get.p, Ytrain)
    pvals_sel = head(pvals[order(pvals$p),],best.dim[i])
    Xtrain1 = Xtrain[,pvals_sel$GENE]
    Xtest1 = Xtest[,pvals_sel$GENE]
    
    model2= naive_bayes(x=Xtrain1,y = Ytrain)
    Ypred = predict(model2,Xtest1)
    cv.error[i]=mean(Ypred!= Ytest)
  }
  tab_naivebayes= cbind(best.dim,cv.error)[1:5,]
  return(tab_naivebayes)
}

##########################################glmnet lasso t-test


lassotst<-function(df,effect)
{
  dim_seq=seq(100,500,50)
  lambda = cv.glmnet(as.matrix(df),as.double(effect))$lambda
  lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
  n <- nrow(expr.1001)
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
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
        
        pvals <- data.frame(GENE = colnames(df))
        pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
        pvals_sel = head(pvals[order(pvals$p),],dim_seq[k])
        #cv.features[j,] = pvals_sel$GENE
        Xtrain2 = Xtrain2[,pvals_sel$GENE]
        Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
        Yvalidate = Ytrain[-train2.ix]
        
        
        fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = lambda_seq[m])
        Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
        inside.error[j] <- mean(Ypred != Yvalidate)
        
      }
      
      lambda.error[k,m] =mean(inside.error)
    }
  }
  # find the k with the smallest error
  idx <- which.min(cost.error)
  if(idx%%9!=0)
  {best_dim=dim_seq[idx%%9]
  best_cost=cost_seq[(idx%/%9+1)]
  }else
  {best_dim=dim_seq[9]
  best_cost=cost_seq[(idx%/%9)]}
  ##best.lambda[i,k] <- lambda_seq[idx]
  pvals <- data.frame(GENE = colnames(df))
  pvals$p = apply(Xtrain,2,get.p, Ytrain)
  pvals_sel = head(pvals[order(pvals$p),],best_dim)
  Xtrain1 = Xtrain[,pvals_sel$GENE]
  Xtest1 = Xtest[,pvals_sel$GENE]
  
  fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda )
  Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
  cv.error <- mean(Ypred!= Ytest)
  tab_glmnet=c(best_lambda,best_dim,cv.error)
  return(tab_glmnet) 
  
}

##########################knn t test

knntst<-function(df,effect)
{
  library(caret)
  
  
  dim_seq=seq(100,500,50)
  lambda_seq=seq(1,19,2)
  # create the vector of tuning parameters
  n <- nrow(expr.1001)
  
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
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
        
        
        pvals <- data.frame(GENE = colnames(expr))
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
      
      lambda.error[i,m] =mean(inside.error)
    }
  }
  # find the k with the smallest error
  idx <- which.min(cost.error)
  if(idx%%9!=0)
  {best_dim=dim_seq[idx%%9]
  best_lambda=lambda_seq[(idx%/%9+1)]
  }else
  {best_dim=dim_seq[9]
  best_lambda=lambda_seq[(idx%/%9)]}
  pvals <- data.frame(GENE = colnames(df))
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
  return(tab_knn)
}


##################################Randomforest

rftst<-function(df,effect)
{
  library(randomForest)
  n=nrow(df)
  dim_seq=seq(500,2000,500)
  dim_seq
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
  dim.error <- vector(mode = "numeric")
  
  for(m in 1:4){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      pvals <- data.frame(GENE = colnames(df))
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
  idx <- which.min(dim.error)
  best.dim<- dim_seq[idx]
  pvals <- data.frame(GENE = colnames(df))
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
  return(tab_rf)
}
 

####################################################svm t test

svmtst<-function(df,effect)
{
  library(e1071)
  
  
  dim_seq=seq(50,500,50)
  cost_seq=seq(0.25,2.5,0.25)
  cost_seq
  # create the vector of tuning parameters
  n <- nrow(df)
  
  
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = df[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = df[-train.ix,]
  Ytest = effect[-train.ix]
  
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
        
        
        pvals <- data.frame(GENE = colnames(df))
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
  if(idx%%10!=0)
  {best_dim=dim_seq[idx%%10]
  best_cost=cost_seq[(idx%/%10+1)]}else
  {best_dim=dim_seq[10]
  best_cost=cost_seq[(idx%/%10)]}
  pvals <- data.frame(GENE = colnames(df))
  pvals$p = apply(Xtrain,2,get.p, Ytrain)
  pvals_sel = head(pvals[order(pvals$p),],best_dim)
  Xtrain1 = Xtrain[,pvals_sel$GENE]
  Xtest1 = Xtest[,pvals_sel$GENE]
  
  
  svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="linear",cost=best_cost)
  Ypred<-predict(svm_model,Xtest1)
  cv.error <- mean(Ypred!= Ytest)
  tab_svm=c(best_cost,best_dim,cv.error)
  return(tab_svm)
        
}

#############################ridge t test
ridge<-function(df,effect)
{dim_seq=seq(100,500,50)
lambda = cv.glmnet(as.matrix(df),as.double(effect),alpha=0)$lambda
lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
n <- nrow(expr.1001)

train.ix <- sample(1:n, (4/5)*n)
Xtrain = df[train.ix,]
Ytrain = effect[train.ix]
Xtest = df[-train.ix,]
Ytest = effect[-train.ix]

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
      
      pvals <- data.frame(GENE = colnames(df))
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
if(idx%%9!=0)
{best_dim=dim_seq[idx%%9]
best_cost=cost_seq[(idx%/%9+1)]
}else
{best_dim=dim_seq[9]
best_cost=cost_seq[(idx%/%9)]}
##best.lambda[i,k] <- lambda_seq[idx]
pvals <- data.frame(GENE = colnames(df))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda ,alpha=0)
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_glmnet=c(best_lambda,best_dim,cv.error)
return(tab_glmnet) 
  
}

##########################################################pcr

pcr<-function(df,effect)             ##############df should be centerd one
{
  library(pls)
  n=dim(df)[1]
  n
  dim_seq=seq(5,n/5,5)
  pvals <- data.frame(GENE = colnames(df))
  pvals$p = apply(center.expr,2,get.p, effect)
  pvals_sel = head(pvals[order(pvals$p),],2000)
  center2000=df[,pvals_sel$GENE]
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = center2000[train.ix,]
  Ytrain = effect[train.ix]
  Xtest = center2000[-train.ix,]
  Ytest = effect[-train.ix]
  train.dat <- data.frame(class = Ytrain, Xtrain)
  test.dat <- data.frame(class =Ytest, Xtest)
  dim.error <- vector(mode = "numeric")
  
  
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
  
  
  tab_pcr=c(best.dim,cv.error)
  
  return(tab_pcr)
  
}


###############################lda t test

ldatst<-function(df,effect)
{
  n=nrow(df)
  cv.error <- vector(mode = "numeric", length = 5)
  best.dim <- vector(mode = "numeric", length = 5)
  dim_seq=seq(50,500,50)
  dim_seq
  for(i in 1:5){
    
    train.ix <- sample(1:n, (4/5)*n)
    Xtrain = df[train.ix,]
    Ytrain = effect[train.ix]
    Xtest = df[-train.ix,]
    Ytest = effect[-train.ix]
    
    dim.error <- vector(mode = "numeric")
    
    for(m in dim_seq){
      
      inside.error <- vector(mode = "numeric", length = 5)
      #cv.features <- matrix(NA,nrow = 5,ncol = 100)
      
      for(j in 1:5){
        n2 <- nrow(Xtrain)
        train2.ix <- sample(1:n2, (4/5)*n2)
        Xtrain2 = Xtrain[train2.ix,]
        Ytrain2 = Ytrain[train2.ix]
        
        pvals <- data.frame(GENE = colnames(df))
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
    pvals <- data.frame(GENE = colnames(df))
    pvals$p = apply(Xtrain,2,get.p, Ytrain)
    pvals_sel = head(pvals[order(pvals$p),],best.dim[i])
    Xtrain1 = Xtrain[,pvals_sel$GENE]
    Xtest1 = Xtest[,pvals_sel$GENE]
    
    model2= lda(Xtrain1,Ytrain)
    Ypred = predict(model2,Xtest1)$class
    cv.error[i]=mean(Ypred!= Ytest)
  }
  tab_lda= cbind(best.dim,cv.error)[1:5,]
  return(tab_lda)
}

############################qda t test
qdatst<-function(df,effect)
{
  n=nrow(df)
  cv.error <- vector(mode = "numeric", length = 5)
  best.dim <- vector(mode = "numeric", length = 5)
  dim_seq=seq(50,500,50)
  dim_seq
  for(i in 1:5){
    
    train.ix <- sample(1:n, (4/5)*n)
    Xtrain = df[train.ix,]
    Ytrain = effect[train.ix]
    Xtest = df[-train.ix,]
    Ytest = effect[-train.ix]
    
    dim.error <- vector(mode = "numeric")
    
    for(m in dim_seq){
      
      inside.error <- vector(mode = "numeric", length = 5)
      #cv.features <- matrix(NA,nrow = 5,ncol = 100)
      
      for(j in 1:5){
        n2 <- nrow(Xtrain)
        train2.ix <- sample(1:n2, (4/5)*n2)
        Xtrain2 = Xtrain[train2.ix,]
        Ytrain2 = Ytrain[train2.ix]
        
        pvals <- data.frame(GENE = colnames(df))
        pvals$p = apply(Xtrain2,2,get.p, Ytrain2)
        pvals_sel = head(pvals[order(pvals$p),],m)
        #cv.features[j,] = pvals_sel$GENE
        Xtrain2 = Xtrain2[,pvals_sel$GENE]
        Xvalidate = Xtrain[-train2.ix,pvals_sel$GENE]
        Yvalidate = Ytrain[-train2.ix]
        
        model1 = qda(Xtrain2,Ytrain2)
        Ypred = predict(model1,Xvalidate)$class
        inside.error[j] <- mean(Ypred != Yvalidate)
      }
      
      dim.error <- append(dim.error, mean(inside.error))
    }
    idx <- which.min(dim.error)
    best.dim[i] <- dim_seq[idx]
    pvals <- data.frame(GENE = colnames(df))
    pvals$p = apply(Xtrain,2,get.p, Ytrain)
    pvals_sel = head(pvals[order(pvals$p),],best.dim[i])
    Xtrain1 = Xtrain[,pvals_sel$GENE]
    Xtest1 = Xtest[,pvals_sel$GENE]
    
    model2= qda(Xtrain1,Ytrain)
    Ypred = predict(model2,Xtest1)$class
    cv.error[i]=mean(Ypred!= Ytest)
  }
  tab_lda= cbind(best.dim,cv.error)[1:5,]
  return(tab_lda)
}


