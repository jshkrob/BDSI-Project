
library(tidyverse)
library(readr)

#PCA IN R
#Once again subset data for drug 1014 only
screen.dat <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/screening.rds")
meth.dat <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/methylation.rds")
effective.1014 <- subset(screen.dat, DRUG_ID_lib == "1014")[,c("CL","CELL_LINE_NAME","EFFECT")]
dim(effective.1014)
meth.1014 <- meth.dat[as.character(effective.1014$CELL_LINE_NAME),]
rm(screen.dat,meth.dat)

#Same t-test function as we used for gene expression data
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

#Ordering pvalues in order of smallest to largest, then selecting the first 10000 p-values.  This is equivalent to the 10000 covariates with maximum variation
pvals <- data.frame(METHSITE = colnames(meth.1014))
pvals$p = apply(meth.1014,2,get.p, effective.1014$EFFECT)

pvals = pvals[complete.cases(pvals),] ### remove values resulting with NAs.
##We will talk of dealing with missing values later.
pvals_sel = head(pvals[order(pvals$p),],10000)

meth.1014.1 = meth.1014[,as.character(pvals_sel$METHSITE)]

#Checkpoint, save the subsetted data as RDS
rm(meth.1014)
saveRDS(meth.1014.1,"./methylation_sel.rds")

meth <- readRDS("./methylation_sel.rds")
dim(meth)

# first, we need to center the expression data by column means
center.meth = scale(meth, scale = F) #scale = F centers the data
# now use svd to perform PCA
pca <- svd(center.meth, nu = 20, nv = 0)
hist(log(pca$d))## plots log of distribution of sigular values (variance of the singular components)

m = sum(cumsum(pca$d^2)/sum(pca$d^2)<.8) #Looks at the cumulative sum vector of the variances over the sum of variance and finds the number of components to explain 80% of the variation
pc.scores <- as.data.frame(diag(pca$d)%*%pca$u) #computes pc score
pc.scores <- cbind(effective.1014, pc.scores) #matches to labels of data
m
dim(pc.scores)

head(pc.scores)

options(repr.plot.width = 5, repr.plot.height = 3)
p1 = ggplot(data = pc.scores) + geom_point(aes(x = V1, y = V2, color = EFFECT)) + labs(x = "PC1", y = "PC2")
p1

p2 = ggplot(data = pc.scores) + geom_point(aes(x = V1, y = V3, color = EFFECT)) + labs(x = "PC1", y = "PC3")
p2

p3 = ggplot(data = pc.scores) + geom_point(aes(x = V3, y = V2, color = EFFECT)) + labs(x = "PC3", y = "PC2")
p3

saveRDS(pc.scores,"./meth_pcscores.rds")

#install.packages("pls")
library(pls)

inside_errors <- vector(mode='numeric', length=10)
set.seed(23527)
for (j in 1:10) {
    n <- dim(center.meth)[1]
    train.index <- sample(1:n, round(.9*n))

    # training set

    x.train <- center.meth[train.index,]
    y.train <- effective.1014[train.index,]
    train.dat <- data.frame(class = y.train$EFFECT, x.train)

    # testing set
    x.test <- center.meth[-train.index,]
    y.test <- effective.1014[-train.index,]
    test.dat <- data.frame(class = y.test$EFFECT, x.test)

    pcr.fit <- pcr(class ~ ., data = train.dat, scale = TRUE, ncomp = 20) #20 principal components so ncomp=20

    # now use that fit to get predictions for the testing data
    pcr.pred <- predict(pcr.fit, test.dat, ncomp = 20)
    preds <- ifelse(pcr.pred >= 0.5, 1, 0)

    # get the classification success rate
    inside_errors[j] <- mean(preds != y.test$EFFECT)    
}
print(mean(inside_errors))


#install.packages('randomForest')
library(randomForest)

# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=10)

n <- dim(center.meth)[1]


# training set
for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
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

print(mean(errors))

#install.packages('caret')
library(caret)

#Installed for CV purposes
# install.packages("class")
library("class")

#Find k which minimizes CV error
set.seed(23527)
M <- seq(1,20,2)    # create the vector of tuning parameters
n <- nrow(pc.scores)

# create a vector to hold the 10 error rates
cv.error <- vector(mode = "numeric", length = 10)

# create a vector to hold the 10 best values of k
best.k <- vector(mode = "numeric", length = 10)

for(i in 1:10){
    
    train.ix <- sample(1:n, (9/10)*n)
    meth.train = pc.scores[train.ix,-c(1:3)]
    effect.train = pc.scores[train.ix,'EFFECT']
    meth.test = pc.scores[-train.ix,-c(1:3)]
    effect.test = pc.scores[-train.ix,'EFFECT']
    
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
    
    fit.train <- knn(train = meth.train, test = meth.test, cl = effect.train, k = M[idx])
    cv.error[i] <- mean(fit.train != effect.test)
    
}
cbind(best.k,cv.error)

# Try k=5
# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=20)

n <- dim(center.meth)[1]

for(j in 1:20){
    train.index <- sample(1:n, round(.95*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
            
    fit <- knn(train = meth.train2, test = meth.validate, cl = effect.train2, k = 5)
    errors[j] <- mean(fit != effect.validate)
}
print(mean(errors))

#Try k=7

# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=20)

n <- dim(center.meth)[1]

for(j in 1:20){
    train.index <- sample(1:n, round(.95*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
            
    fit <- knn(train = meth.train2, test = meth.validate, cl = effect.train2, k = 7)
    errors[j] <- mean(fit != effect.validate)
}
print(mean(errors))

#install.packages("e1071")
library(e1071)


# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=10)

n <- nrow(pc.scores)


# training set
for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    #dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
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

print(mean(errors))

# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=10)

n <- nrow(pc.scores)


# training set
for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    #dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
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

print(mean(errors))

# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=10)

n <- nrow(pc.scores)


# training set
for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    #dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
    test.dat <- data.frame(class = y.test, x.test)
    #dim(test.dat)
    svm_model <- svm(as.factor(y.train) ~ ., data=x.train, kernel="sigmoid")
    #Ytrainhat = rf.fit$predicted
    preds <- predict(svm_model,x.test)
    #trainingerror$rf = 1 - mean(Ytrainhat==Ytrain)
    #testerror$rf = 1 - mean(Ytesthat==Ytest)

    # get the classification success rate
    errors[i] <- mean(preds != y.test)    
}

print(mean(errors))

#install.packages('HiDimDA')
library(HiDimDA)

# split the centered data into training and testing
set.seed(23527)
errors <- vector(mode='numeric', length=10)

n <- nrow(pc.scores)


# training set
for (i in 1:10) {
    train.index <- sample(1:n, round(.9*n))
    x.train <- pc.scores[train.index,-c(1,2,3)]
    y.train <- pc.scores[train.index,'EFFECT']
    train.dat <- data.frame(class = y.train, x.train)
    #dim(train.dat)
    # testing set
    x.test <- pc.scores[-train.index,-c(1,2,3)]
    y.test <- pc.scores[-train.index,'EFFECT']
    test.dat <- data.frame(class = y.test, x.test)
    #dim(test.dat)
    dlda_model <- Dlda(data=train.dat, grouping=as.factor(train.dat$class), VSelfunct = "none",
        ldafun="classification")
    #Ytrainhat = rf.fit$predicted
    preds <- predict(dlda_model,test.dat)$class
    #print(summary(preds))
    #trainingerror$rf = 1 - mean(Ytrainhat==Ytrain)
    #testerror$rf = 1 - mean(Ytesthat==Ytest)

    # get the classification success rate
    errors[i] <- mean(preds != y.test)    
}

print(mean(errors))
