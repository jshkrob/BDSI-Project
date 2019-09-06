screening<-readRDS("screening.rds")
expr=readRDS("expression.rds")
effective.1066<- subset(screening, DRUG_ID_lib == "1066")[,c("CL","EFFECT")]
expr.1066 <- expr[as.character(effective.1066$CL),]
Y=effective.1066$EFFECT
Y
center.expr=scale(expr.1066,scale=F)
library(matlib)


meth=readRDS("methylation_na.rds")
screening<-readRDS("screening.rds")
effective.1066 <- subset(screening, DRUG_ID_lib == "1066")[,c("CL","EFFECT","CELL_LINE_NAME")]
mds.1066 = meth[as.character(effective.1066$CELL_LINE_NAME),]
m=20

pccomp=prcomp(center.expr)

pcexpr=pccomp$x[,1:m]
a=pcexpr
proj=a%*%Ginv(t(a)%*%a)%*%t(a)

orthoproj=(diag(nrow(expr.1066))-proj)
orthometh=orthoproj%*%as.matrix(mds.1066)
dim(orthometh)

center.meth=scale(orthometh,scale=F)
pccomp1=prcomp(center.meth)  
pcmeth=pccomp1$x[,1:m]
b=pcmeth

c=cbind(a,b)

proj1=c%*%Ginv(t(c)%*%c)%*%t(c)

dim(proj1)
orthoproj1=(diag(nrow(expr.1066))-proj1)


R=orthoproj1

#####################################################################################
copynum=readRDS("copynumber.rds")
sum(is.na(copynum))
dim(copynum)
copynum=copynum[,colSums(is.na(copynum))<850]
copy.1066 <- copynum[as.character(effective.1066$CL),]
n=nrow(copy.1066)
scale=matrix(0,nrow=n,ncol=n)
C=copy.1066
for(i in 1:n)
{for(j in 1:n)
{k=sum(is.na(C[i,]+C[j,]))
scale[i,j]=ncol(C)/(ncol(C)-k)
}
}

C[is.na(C)]<-0
sum(is.na(C))

C1<-C%*%t(C)
dim(C1)

C1<-C1*scale
C1

C2<-R%*%C1%*%R
C2<-scale(C2,scale=F)

sv=svd(C2,nu=m,nv=0)

d=diag(sv$d[1:m])
pccopy=sv$u%*%sqrt(d)
data=cbind(pcexpr,pcmeth,pccopy)



dim(data)

colnames(data)=as.character(c(1:(3*m)))
m

########################Lda
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



###################################qda
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
########################################NaiveBAyes

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

###############################randomForest

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

##################################logistic glmnet lasso

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



#######################glmnet ridge

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


############################svm

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
###########################knn

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
combqda(data,Y)
comblda(data,Y)
combknn(data,Y)
combsvm(data,Y)
combridge(data,Y)
comblasso(data,Y)
combrf(data,Y)
combnb(data,Y)

