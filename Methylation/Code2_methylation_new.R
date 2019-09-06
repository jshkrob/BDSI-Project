meth=readRDS("methylation_na.rds")
sum(is.na(meth))
dim(meth)
effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT","CELL_LINE_NAME")]
mds.1014 = meth[as.character(effective.1014$CELL_LINE_NAME),]
dim(mds.1014)
center.meth=scale(mds.1014)
n=nrow(effective.1014)
pca_comp=prcomp(center.meth)
pc.scoresdat<-pca_comp$x
dim_seq=seq(5,60,5)
dim_seq
sum(is.na(pc.scoresdat))
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)
dim(pc.scoresdat)
best.dim

library(randomForest)
for(i in 1:5)
{ 
train.ix <- sample(1:n, 0.8*n)
  Xtrain = pc.scoresdat[train.ix,]
  Ytrain = effective.1014$EFFECT[train.ix]
  Xtest = pc.scoresdat[-train.ix,]
  Ytest = effective.1014$EFFECT[-train.ix]
  train = data.frame(Xtrain, Ytrain)      # assemble dataframes
  test = data.frame(Xtest, Ytest)
  dim.error<-vector(mode="numeric")
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,1:m]
      Ytrain2 = Ytrain[train2.ix]
      Xvalidate = Xtrain[-train2.ix,1:m]
      Yvalidate = Ytrain[-train2.ix]
      train2 = data.frame(Xtrain2, Ytrain2)      # assemble dataframes
      test2 = data.frame(Xvalidate, Yvalidate)
      
      model6 <- randomForest(as.factor(Ytrain2) ~ ., data=train2)
      Ypred = predict(model6,Xvalidate)
      inside.error[j] <- mean(Ypred != Yvalidate)
    }
    
    dim.error <- append(dim.error, mean(inside.error))
  }
  idx <- which.min(dim.error)
  best.dim[i]<- dim_seq[idx]
  idx
  best.dim[1]
  Xtrain1 = Xtrain[,1:best.dim[i]]
  Xtest1 = Xtest[,1:best.dim[i]]
  train1=data.frame(Xtrain1,Ytrain)
  model6 <- randomForest(as.factor(Ytrain) ~ ., data=train1)
  Ypred = predict(model6,Xtest1)
  cv.error[i]=mean(Ypred!= Ytest)
}
tab_pcrf= cbind(best.dim,cv.error)[1:5,]
tab_pcrf

##############pcalasso
library(glmnet)


dim_seq=seq(5,60,5)
dim_seq
lambda = cv.glmnet(as.matrix(pc.scoresdat),as.double(effective.1014$EFFECT))$lambda

lambda_seq = seq(min(lambda),max(lambda),length.out = 10) 
lambda_seq
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

lambda.error <- matrix(0,nrow=12,ncol=10)
lambda.error
for(i in 1:12){
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
      
      fit <- glmnet(x = as.matrix(Xtrain2), y = as.double(Ytrain2),family = "binomial",lambda = lambda_seq[m],alpha=1)
      Ypred = predict(fit,as.matrix(Xvalidate),type = "response")>.5
      inside.error[j] <- mean(Ypred != Yvalidate)
      
    }
    
    lambda.error[i,m] =mean(inside.error)
  }
}
# find the k with the smallest error
idx <- which.min(lambda.error)
if(idx%%12!=0)
{best_dim=dim_seq[idx%%12]
best_lambda=lambda_seq[(idx%/%12+1)]
}else
{best_dim=dim_seq[12]
best_lambda=lambda_seq[(idx%/%12)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda )
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_lasso=c(best_lambda,best_dim,cv.error)
t

######################svm

library(e1071)


dim_seq=seq(5,60,5)
cost_seq=seq(0.25,2.5,0.25)
# create the vector of tuning parameters
n <- nrow(pc.scoresdat)
n


train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

cost.error <- matrix(0,nrow=12,ncol=10)
cost.error
for(i in 1:12){
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
if(idx%%12!=0)
{best_dim=dim_seq[idx%%12]
best_cost=cost_seq[(idx%/%12+1)]
}else
{best_dim=dim_seq[12]
best_cost=cost_seq[(idx%/%12)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="linear",cost=best_cost)
Ypred = predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svm=c(best_cost,best_dim,cv.error)
tab_svm


###############knn

library(caret)


dim_seq=seq(5,60,5)
lambda_seq=seq(1,19,2)
# create the vector of tuning parameters
n <- nrow(pc.scoresdat)



train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

lambda.error <- matrix(0,nrow=12,ncol=10)
lambda.error
for(i in 1:12){
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
if(idx%%12!=0)
{best_dim=dim_seq[idx%%12]
best_lambda=lambda_seq[(idx%/%12+1)]
}else
{best_dim=dim_seq[12]
best_lambda=lambda_seq[(idx%/%12)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

model5a = knn3(Ytrain ~ ., data=train1,  k = best_lambda) 
Ypred = ifelse(predict(model5a,data.frame(Xtest1))[,2]>=.5,1,0)
cv.error <- mean(Ypred!= Ytest)
tab_knn=c(best_lambda,best_dim,cv.error)
tab_knn


################naivebayes
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)

dim_seq=seq(5,60,5)
library(naivebayes)
for(i in 1:5)
{train.ix <- sample(1:n, 0.8*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]

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


######################lda

cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)

dim_seq=seq(5,60,5)
library(MASS)
for(i in 1:5)
{train.ix <- sample(1:n, 0.8*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]

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


################################qda
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)

dim_seq=seq(5,60,5)
library(MASS)
for(i in 1:5)
{train.ix <- sample(1:n, 0.8*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]

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
    
    model1 = qda(Xtrain2,Ytrain2)
    Ypred = predict(model1,Xvalidate)$class
    inside.error[j] <- mean(Ypred != Yvalidate)
  }
  
  dim.error <- append(dim.error, mean(inside.error))
}
idx <- which.min(dim.error)
best.dim[i]<- dim_seq[idx]

Xtrain1 = Xtrain[,1:best.dim[i]]
Xtest1 = Xtest[,1:best.dim[i]]

model2= qda(Xtrain1,Ytrain)
Ypred = predict(model2,Xtest1)$class
cv.error[i]=mean(Ypred!= Ytest)
}
tab_pcqda= cbind(best.dim,cv.error)[1:5,]
tab_pcqda

################################glmnetridge

library(glmnet)


dim_seq=seq(5,60,5)
lambda = cv.glmnet(as.matrix(pc.scoresdat),as.double(effective.1014$EFFECT),alpha=0)$lambda

lambda_seq = seq(min(lambda),max(lambda),length.out = 10) 
lambda_seq
train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]
train = data.frame(Xtrain, Ytrain)      # assemble dataframes
test = data.frame(Xtest, Ytest)

lambda.error <- matrix(0,nrow=12,ncol=10)
lambda.error
for(i in 1:12){
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
if(idx%%12!=0)
{best_dim=dim_seq[idx%%12]
best_lambda=lambda_seq[(idx%/%12+1)]
}else
{best_dim=dim_seq[12]
best_lambda=lambda_seq[(idx%/%12)]}
Xtrain1 = Xtrain[,1:best_dim]
Xtest1 = Xtest[,1:best_dim]
train1 = data.frame(Xtrain1, Ytrain)      # assemble dataframes
test1 = data.frame(Xtest1, Ytest)

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda ,alpha=0)
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_ridge=c(best_lambda,best_dim,cv.error)
tab_ridge


#########################pcr
library(pls)
n=dim(center.meth)[1]
n
dim_seq=seq(5,60,5)
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

pvals <- data.frame(GENE = colnames(center.meth))
pvals$p = apply(center.meth,2,get.p, Y)
dim(center.meth)
pvals_sel = head(pvals[order(pvals$p),],2000)
center.meth1000=center.meth[,pvals_sel$GENE]
train.ix <- sample(1:n, (4/5)*n)
Xtrain = center.meth1000[train.ix,]
Ytrain = effective.1014$EFFECT[train.ix]
Xtest = center.meth1000[-train.ix,]
Ytest = effective.1014$EFFECT[-train.ix]
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


tab_pcr=cbind(best.dim,cv.error)

tab_pcr
