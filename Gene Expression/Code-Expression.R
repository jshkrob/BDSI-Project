screening<-readRDS("screening.rds")
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


#############################DRUG1001

expr=readRDS("expression.rds")
effective.1001 <- subset(screening, DRUG_ID_lib == "1001")[,c("CL","EFFECT")]
expr.1001 <- expr[as.character(effective.1001$CL),]
Y=effective.1001$EFFECT
Y
dim(effective.1001)
dim(expr.1001)

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
################naivebayes
n=nrow(expr.1001)
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)
dim_seq=seq(50,500,50)
dim_seq
for(i in 1:5){
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = expr.1001[train.ix,]
  Ytrain = effective.1001$EFFECT[train.ix]
  Xtest = expr.1001[-train.ix,]
  Ytest = effective.1001$EFFECT[-train.ix]
  
  dim.error <- vector(mode = "numeric")
  
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
  pvals <- data.frame(GENE = colnames(expr.1001))
  pvals$p = apply(Xtrain,2,get.p, Ytrain)
  pvals_sel = head(pvals[order(pvals$p),],best.dim[i])
  Xtrain1 = Xtrain[,pvals_sel$GENE]
  Xtest1 = Xtest[,pvals_sel$GENE]
  
  model2= naive_bayes(x=Xtrain1,y = Ytrain)
  Ypred = predict(model2,Xtest1)
  cv.error[i]=mean(Ypred!= Ytest)
}
tab_naivebayes= cbind(best.dim,cv.error)[1:5,]
tab_naivebayes
testerror=c(rep(0,times=10))
testerror[1]=(cv.error[1]+cv.error[4])/2
testerror
bestdim_nb=350

######################################################lda
n=nrow(expr.1001)
##cv.error <- vector(mode = "numeric", length = 5)
##best.dim <- vector(mode = "numeric", length = 5)
dim_seq=seq(50,500,50)
dim_seq
for(i in 1:5){
  
  train.ix <- sample(1:n, (4/5)*n)
  Xtrain = expr.1001[train.ix,]
  Ytrain = effective.1001$EFFECT[train.ix]
  Xtest = expr.1001[-train.ix,]
  Ytest = effective.1001$EFFECT[-train.ix]
  
  dim.error <- vector(mode = "numeric")
  
  for(m in dim_seq){
    
    inside.error <- vector(mode = "numeric", length = 5)
    #cv.features <- matrix(NA,nrow = 5,ncol = 100)
    
    for(j in 1:5){
      n2 <- nrow(Xtrain)
      train2.ix <- sample(1:n2, (4/5)*n2)
      Xtrain2 = Xtrain[train2.ix,]
      Ytrain2 = Ytrain[train2.ix]
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
  pvals <- data.frame(GENE = colnames(expr.1001))
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
bestdim_lda=450


################################glmnet
library(glmnet)


dim_seq=seq(100,500,50)
lambda = cv.glmnet(as.matrix(expr.1001),as.double(effective.1001$EFFECT))$lambda
lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
n <- nrow(expr.1001)


# create a vector to hold the 10 error rate

# create a vector to hold the 10 best values of k
##best.lambda <-matrix(0,nrow=9,ncol=10)
##finerror= vector(mode = "numeric", length = 5)
##finlambda= vector(mode = "numeric", length = 5)
##findim= vector(mode = "numeric", length = 5)
#best.features <- matrix(NA,nrow = 10,ncol = 100)


###for(i in 1:5){
train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]

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
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
pvals <- data.frame(GENE = colnames(expr.1001))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda )
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_glmnet=c(best_lambda,best_dim,cv.error)
tab_glmnet  

##finerror[i]=min(cv.error[i,])
##finlambda[i]=best.lambda(which.min(cv.error[i,]))
##findim[i]=(which.min(cv.error[i,]))*100
##}

##tab=cbind(findim,finerror,finlambda)[1:5,]


###################################knn
library(caret)


dim_seq=seq(100,500,50)
lambda_seq=seq(1,19,2)
# create the vector of tuning parameters
n <- nrow(expr.1001)



train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]
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
      
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
pvals <- data.frame(GENE = colnames(expr.1001))
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
## n=nrow(expr.1001)
##library(randomForest)
##pvals <- data.frame(GENE = colnames(expr.1001))
##pvals$p = apply(expr.1001,2,get.p, Y)
##pvals_sel = head(pvals[order(pvals$p),],10000)


#inside.error <- vector(mode = "numeric", length = 5)
##for(i in 1:5)
##{train.ix <- sample(1:n, (4/5)*n)
# Xtrain = expr.1001[train.ix,]
#Ytrain = effective.1001$EFFECT[train.ix]
#Xtest = expr.1001[-train.ix,]
#Ytest = effective.1001$EFFECT[-train.ix]
#Xtrain1 = Xtrain[,pvals_sel$GENE]
#Xtest1 = Xtest[,pvals_sel$GENE]
#train1 = data.frame(Xtrain1, Ytrain)
#model6 <- randomForest(as.factor(Ytrain) ~ ., data=train1)
#Ypred = predict(model6,Xtest1)
#inside.error[i]=mean(Ypred!= Ytest)

#}

#cv.error=mean(inside.error)
#cv.error 

n=nrow(expr.1001)
dim_seq=seq(500,5000,500)
dim_seq


train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]

dim.error <- vector(mode = "numeric")

for(m in 1:10){
  
  inside.error <- vector(mode = "numeric", length = 5)
  #cv.features <- matrix(NA,nrow = 5,ncol = 100)
  
  for(j in 1:5){
    n2 <- nrow(Xtrain)
    train2.ix <- sample(1:n2, (4/5)*n2)
    Xtrain2 = Xtrain[train2.ix,]
    Ytrain2 = Ytrain[train2.ix]
    
    pvals <- data.frame(GENE = colnames(expr.1001))
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
pvals <- data.frame(GENE = colnames(expr.1001))
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
###########################Svm
library(e1071)


dim_seq=seq(50,500,50)
cost_seq=seq(0.25,2.5,0.25)
cost_seq
# create the vector of tuning parameters
n <- nrow(expr.1001)



train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]

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
      
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
pvals <- data.frame(GENE = colnames(expr.1001))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]


svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="linear",cost=best_cost)
Ypred<-predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svm=c(best_cost,best_dim,cv.error)
tab_svm
warnings()


#################################################glmnet-ridge

library(glmnet)


dim_seq=seq(100,500,50)
lambda = cv.glmnet(as.matrix(expr.1001),as.double(effective.1001$EFFECT),alpha=0)$lambda
lambda_seq = seq(min(lambda),max(lambda),length.out = 10)    # create the vector of tuning parameters
n <- nrow(expr.1001)
lambda



# create a vector to hold the 10 error rate

# create a vector to hold the 10 best values of k
##best.lambda <-matrix(0,nrow=9,ncol=10)
##finerror= vector(mode = "numeric", length = 5)
##finlambda= vector(mode = "numeric", length = 5)
##findim= vector(mode = "numeric", length = 5)
#best.features <- matrix(NA,nrow = 10,ncol = 100)


###for(i in 1:5){
train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]

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
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
best_cost=cost_seq[(idx%/%9+1)]}else
{best_dim=dim_seq[9]
best_cost=cost_seq[(idx%/%9)]}
##best.lambda[i,k] <- lambda_seq[idx]
pvals <- data.frame(GENE = colnames(expr.1001))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]

fit.train <- glmnet(x = as.matrix(Xtrain1), y = as.double(Ytrain),family = "binomial",lambda = best_lambda ,alpha=0)
Ypred = predict(fit.train,as.matrix(Xtest1),type = "response")>.5
cv.error <- mean(Ypred!= Ytest)
tab_glmnetridge=c(best_lambda,best_dim,cv.error)
tab_glmnetridge



##################################################svm-polynomial 

library(e1071)


dim_seq=seq(50,500,50)
cost_seq=seq(0.25,2.5,0.25)
cost_seq
# create the vector of tuning parameters
n <- nrow(expr.1001)



train.ix <- sample(1:n, (4/5)*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]

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
      
      
      pvals <- data.frame(GENE = colnames(expr.1001))
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
if(idx%%10!=0)
{best_dim=dim_seq[idx%%10]
best_cost=cost_seq[(idx%/%10+1)]}else
{best_dim=dim_seq[10]
best_cost=cost_seq[(idx%/%10)]}
pvals <- data.frame(GENE = colnames(expr.1001))
pvals$p = apply(Xtrain,2,get.p, Ytrain)
pvals_sel = head(pvals[order(pvals$p),],best_dim)
Xtrain1 = Xtrain[,pvals_sel$GENE]
Xtest1 = Xtest[,pvals_sel$GENE]


svm_model <- svm(as.factor(Ytrain)~ ., data=Xtrain1, kernel="polynomial",cost=best_cost)
Ypred<-predict(svm_model,Xtest1)
cv.error <- mean(Ypred!= Ytest)
tab_svmpoly=c(best_cost,best_dim,cv.error)
tab_svmpoly
warnings()


##################################################################rda

n=nrow(expr.1001)
library(rda)

train.ix <- sample(1:n, 0.6*n)
Xtrain = expr.1001[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = expr.1001[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]
train=data.frame(Xtrain,Ytrain)

fit<-rda(x = t(as.matrix(Xtrain)),y = as.double(Ytrain), xnew = t(Xtest),ynew = as.double(Ytest),delta = 0.0001,alpha = 0)

# Ypred<-predict(fit,t(Xtrain),Ytrain,t(Xtest))
# error<-mean(Ypred!=Ytest)
# error
##################################################
##################################################PCA
expr=readRDS("expression.rds")
screening = readRDS("screening.rds")
effective.1001 <- subset(screening, DRUG_ID_lib == "1001")[,c("CL","CELL_LINE_NAME","EFFECT")]
dim(effective.1001)
expr.1001 <- expr[as.character(effective.1001$CL),]
center.expr=scale(expr.1001,scale=F)
View(center.expr)
pca <- svd(center.expr, nu = 20, nv = 0)  
pca  
m = sum(cumsum(pca$d^2)/sum(pca$d^2)<0.9)
m
pca_comp=prcomp(center.expr)
pc.scoresdat<-pca_comp$x
dim(pc.scoresdat)
pc.scores <- cbind(effective.1001, pc.scoresdat)
dim(pc.scores)
saveRDS(pc.scores,"pcaexpression-90%variance.rds")

########################################################pcr

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
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = center.expr1000[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]
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



############################pcalda

n=nrow(expr.1001)
pca_comp=prcomp(center.expr)
pc.scoresdat<-pca_comp$x
cv.error <- vector(mode = "numeric", length = 5)
best.dim <- vector(mode = "numeric", length = 5)
train.ix <- sample(1:n, 0.8*n)
dim_seq=seq(5,40,5)
library(MASS)
for(i in 1:5)
{
  Xtrain = pc.scoresdat[train.ix,]
  Ytrain = effective.1001$EFFECT[train.ix]
  Xtest = pc.scoresdat[-train.ix,]
  Ytest = effective.1001$EFFECT[-train.ix]
  
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


#########################################pcaknn
library(caret)


dim_seq=seq(5,40,5)
lambda_seq=seq(1,19,2)
# create the vector of tuning parameters
n <- nrow(expr.1001)



train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]
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

model5a = knn3(Ytrain ~ ., data=train1,  k = best_lambda) 
Ypred = ifelse(predict(model5a,data.frame(Xtest1))[,2]>=.5,1,0)
cv.error <- mean(Ypred!= Ytest)
tab_knn=c(best_lambda,best_dim,cv.error)
tab_knn


################################svmlinearpca

library(e1071)


dim_seq=seq(5,40,5)
cost_seq=seq(0.25,2.5,0.25)
# create the vector of tuning parameters
n <- nrow(pc.scoresdat)



train.ix <- sample(1:n, (4/5)*n)
Xtrain = pc.scoresdat[train.ix,]
Ytrain = effective.1001$EFFECT[train.ix]
Xtest = pc.scoresdat[-train.ix,]
Ytest = effective.1001$EFFECT[-train.ix]
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


