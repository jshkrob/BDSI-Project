  ###################################
#JACOBS STUFF GOES HERE #############
##################################
set.seed(23527)
#FUNCTION FORMAT
#pcasvm_ensemble<-function(df, effect, roc = F, stack = FALSE, train.ix, best.cost, best.dim) 
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
meth <- readRDS("methylation_na.rds")
library(readxl)
# Read in data that going to use for adjusting due to tissue site
cell_line_type <- read_excel("~/Cell_Lines_Details.xlsx", sheet = "Cell line details")
colnames(cell_line_type) <- c("Sample_Name", "CL", "WES", "CNA", "GeneExp", "Methyl", "Drug_Resp",
                              "Site", "Tissue_desc_2", "Cancer_Type", "MSI", "Medium", "GrowthProp")
screening <- screening %>% left_join(cell_line_type, by = "CL")

effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","EFFECT", "CELL_LINE_NAME","Site")]


###-------------------------------------------------------------for specific drug, load expr, meth, copynum data ------------------------------------------------

# Make scale function
scale_this <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  return(x - mu)
}

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

## drug id: 1014


### Expression

expr.1014 <- expr[as.character(effective.1014$CL), ]

expr.1014 <- as.data.frame(expr.1014)

expr.1014 <- expr.1014 %>% 
  add_column(Site = effective.1014$Site)

tissuecenter.expr.1014 <- expr.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% # This removes site so it can run the scale_this function on all other vars
  ungroup() %>% 
  dplyr::select(-Site)

effective.c.1014 <- effective.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  ungroup() %>% 
  dplyr::select(-Site) %>% 
  pull(EFFECT)
  
pvals <- data.frame(GENE = colnames(tissuecenter.expr.1014))
pvals$p = apply(tissuecenter.expr.1014,2,get.p,effective.c.1014)
pvals_sel = head(pvals[order(pvals$p),],10000)
tissuecenter.expr.1014.1 = tissuecenter.expr.1014[,as.character(pvals_sel$GENE)]
center.expr.1014 = scale(tissuecenter.expr.1014, scale=F)
center.expr.1014.1 = scale(tissuecenter.expr.1014.1, scale=F)

### Methylation

mds.1014 <- meth[as.character(effective.1014$CELL_LINE_NAME), ]
mds.1014 <- mds.1014[ , !is.na(colSums(meth))]

mds.1014 <- as.data.frame(mds.1014)
mds.idx <- sample(1:ncol(mds.1014),size=20000) #take a random sample of 20000 columns

mds.1014 <- mds.1014[,mds.idx] %>% 
  add_column(Site = effective.1014$Site)

tissuecenter.mds.1014 <- mds.1014 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% # This removes site so it can run the scale_this function on all other vars
  ungroup() %>% 
  dplyr::select(-Site)


pvals <- data.frame(SITE = colnames(tissuecenter.mds.1014))
pvals$p = apply(tissuecenter.mds.1014,2,get.p,effective.c.1014)
pvals_sel = head(pvals[order(pvals$p),],10000)
tissuecenter.mds.1014.1 = tissuecenter.mds.1014[,as.character(pvals_sel$SITE)]
center.meth.1014 = scale(tissuecenter.mds.1014, scale =F)
center.meth.1014.1 = scale(tissuecenter.mds.1014.1, scale=F)

### Copy Number
copynum.1014 <- readRDS(paste0("~/Mcopynum", 1014, "_0.1_KEEP_RMNA.impDATA.rds"))

# Read in datasets corresponding to specific drug
cn.1014 <- paste0("~/Mcopynum", 1014, "_0.1_KEEP_RMNA.impDATA.rds") %>% readRDS(.)

# Make sure both matrices have same cell lines
CN <- cn.1014[as.character(effective.1014$CL),]

CN <- as.data.frame(CN)

CN <- CN %>% 
  add_column(Site = effective.1014$Site)

#Center by tissue site
tissuecenter.cn.1014 <- CN %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% # This removes site so it can run the scale_this function on all other vars
  ungroup() %>% 
  dplyr::select(-Site)


pvals <- data.frame(cnGENE = colnames(tissuecenter.cn.1014))
pvals$p = apply(tissuecenter.cn.1014,2,get.p, effective.c.1014)

pvals = pvals[complete.cases(pvals),] ### remove values resulting with NAs.
##We will talk of dealing with missing values later.
pvals_sel = head(pvals[order(pvals$p),],10000)

cn.1014.1 = tissuecenter.cn.1014[,as.character(pvals_sel$cnGENE)]

center.cn.1014 = scale(tissuecenter.cn.1014, scale = F) #scale = F centers the data
center.cn.1014.1 = scale(cn.1014.1, scale = F) #scale = F centers the data

#IMPORTANT VARIABLES : 
#center.expr.1014 - expression data without t-test stuff
#center.expr.1014.1 - expression data with t-test
#center.meth.1014 - meth data w/o t-test
#center.meth.1014.1 - meth data w t-test
#center.cn.1014 - CN data w/o t-test
#center.cn.1014.1 - CN data w/ t-test
#effective.c.1014 - effective columns relating to 1014
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
f1 <- function(predicted, expected) {
  predict <- ifelse(predicted >= 0.5, 1, 0)
  total_pos_predicted <- sum(predict)
  
  # Precision: fraction of positive predicted results that are actually positive
  precision <- sum(predict & expected) / total_pos_predicted
  
  # Recall: fraction of actual positive results that are predicted pos
  recall <- sum(predict & expected) / sum(expected)
  f1 <- (2*precision*recall) / (precision + recall)
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


f1 <- function(predicted, expected) {
  predict <- ifelse(predicted >= 0.5, 1, 0)
  total_pos_predicted <- sum(predict)
  
  # Precision: fraction of positive predicted results that are actually positive
  precision <- sum(predict & expected) / total_pos_predicted
  
  # Recall: fraction of actual positive results that are predicted pos
  recall <- sum(predict & expected) / sum(expected)
  f1 <- (2*precision*recall) / (precision + recall)
}


#helper functions
sign_chg <- function(x) {
  ifelse(x == 0, -1, 1)
}
is.positive <- function(x) {
  x > 0
}
# ensemble vote method: Chris
ensemble_vote <- function(effect, roc = F, stack = FALSE, train.ix, meth.fit, expr.fit, cn.fit) {
    Ytest <- effect[-train.ix]
    predic_ensemble <- data.frame(ifelse(meth.fit$predicted == TRUE, 1, 0),
                                ifelse(expr.fit$predicted == TRUE, 1, 0),
                                ifelse(copy.fit$predicted == TRUE, 1, 0))
    colnames(predic_ensemble) <- c("METH", "EXPR", "COPY")
    
    Ypred_ensemble <- predic_ensemble %>% 
    mutate_all(funs(sign_chg(.))) %>% 
    mutate(vote = is.positive(METH + EXPR + COPY)) %>% 
    pull(vote)
    
    ensemble_error <- mean(Ypred_ensemble != Ytest)
    
    vote_roc <- roc(Ytest, as.numeric(Ypred_ensemble))
    
    if (roc==TRUE) {
    
    if (stack == TRUE) {
      plot(vote_roc, col = "#d53e4f", add = TRUE)
      
    } else {
      plot(vote_roc, col = "#d53e4f")
      
    }
    }
    
    AUC <- vote_roc$auc
    F1_SCORE <- f1(Ypred_ensemble, Ytest)
    cv.error=ensemble_error
    vote_sensitivities <- vote_roc$sensitivies
    vote_specificities <- vote_roc$specificities
    vote_specificities <-  1 - vote_specificities
    return(list(results = cv.error, predicted = Ypred_ensemble, auc = AUC, f1score = F1_SCORE, sensitivies = vote_sensitivities, specificities=vote_specificities))
}

### ACCUMULATOR MATRIX
ensemble_classifiers <- matrix(nrow=4,ncol=3)

rownames(ensemble_classifiers) = c('meth','expr','cn','ensemble vote')
colnames(ensemble_classifiers) = c('Avg error','Avg AUC', 'Avg F1')



n = nrow(center.meth.1014.1)
train.ix <- sample(1:n, 0.8*n)

meth.fit <- pcalda_ensemble(center.meth.1014.1,effective.c.1014, roc=F, stack=T,train.ix=train.ix, best.dim=25)
expr.fit <- pcalasso_ensemble(center.expr.1014, effective.c.1014,roc=F,stack=T,train.ix=train.ix,best.lambda=0.003485931,best.dim=45)
cn.fit <- pcalasso_ensemble(center.cn.1014, effective.c.1014,roc=F,stack=T,train.ix=train.ix,best.lambda=0.05814357,best.dim=5)
ensemble.fit <- ensemble_vote(effect = effective.c.1014, roc=F,stack=T, meth.fit=meth.fit, expr.fit=expr.fit, cn.fit=cn.fit, train.ix=train.ix)


#Matrices of format rows are iterations and columns are data types (expr,meth,cn,ensemble)
error <- matrix(0,nrow=5,ncol=4)
colnames(error) = c("expr","meth",'cn','ensemble')

f1 <- matrix(0,nrow=5,ncol=4)
colnames(f1) = c("expr","meth",'cn','ensemble')

auc <- matrix(0,nrow=5,ncol=4)
colnames(auc) = c("expr","meth",'cn','ensemble')

meth_roc <- matrix(c(meth.fit$sensitivities, meth.fit$specificities), nrow=length(meth.fit$sensitivities),ncol=2, byrow=F)
colnames(meth_roc) = c("Sens", "Spec")
expr_roc <- matrix(c(expr.fit$sensitivities, expr.fit$specificities), nrow=length(expr.fit$sensitivities),ncol=2,byrow=F)
colnames(expr_roc) = c("Sens", "Spec")
cn_roc <- matrix(c(cn.fit$sensitivities, cn.fit$specificities), nrow=length(cn.fit$sensitivities),ncol=2,byrow=F)
colnames(cn_roc) = c("Sens", "Spec")
ensemble_roc <- matrix(c(ensemble.fit$sensitivities, ensemble.fit$specificities), nrow=length(ensemble.fit$sensitivities),ncol=2,byrow=F)
colnames(ensemble_roc) = c("Sens", "Spec")

error[1,1] <- expr.fit$results[1]
error[1,2] <- meth.fit$results[1]
error[1,3] <- cn.fit$results[1]
error[1,4] <- ensemble.fit$results

f1[1,1] <- expr.fit$f1score
f1[1,2] <- meth.fit$f1score
f1[1,3] <- cn.fit$f1score
f1[1,4] <- ensemble.fit$f1score

auc[1,1] <- expr.fit$auc
auc[1,2] <- meth.fit$auc
auc[1,3] <- cn.fit$auc
auc[1,4] <- ensemble.fit$auc



for (i in 2:5) { #iterates from 2 to 5 because the first of five iterations is done outside the for loop
  train.ix <- sample(1:n, 0.8*n)
  
  #switch out lda with the classifier from the excel sheet (xl_param and xl_dim will be the same)
  meth.fit <- pcalda_ensemble(center.meth.1014.1,effective.c.1014, roc=F, stack=T,train.ix=train.ix, best.dim=25)
  expr.fit <- pcalasso_ensemble(center.expr.1014, effective.c.1014,roc=F,stack=T,train.ix=train.ix,best.lambda=0.003485931,best.dim=expr.xl_dim)
  cn.fit <- pcalasso_ensemble(center.cn.1014, effective.c.1014,roc=F,stack=F,train.ix=train.ix,best.lambda=0.05814357,best.dim=5)
  ensemble.fit <- ensemble_vote(effect = effective.c.1014, roc=F,stack=T, meth.fit=meth.fit, expr.fit=expr.fit, cn.fit=cn.fit, train.ix=train.ix)
  
  meth_roc[,1] <- meth_roc[,1] + meth.fit$sensitivities
  expr_roc[,1] <- expr_roc[,1] + expr.fit$sensitivities
  cn_roc[,1] <- cn_roc[,1] + cn.fit$sensitivities
  ensemble_roc[,1] <- ensemble_roc[,1] + ensemble.fit$sensitivities
  
  meth_roc[,2] <- meth_roc[,2] + meth.fit$specificities
  expr_roc[,2] <- expr_roc[,2] + expr.fit$specificities
  cn_roc[,2] <- cn_roc[,2] + cn.fit$specificities
  ensemble_roc[,2] <- ensemble_roc[,2] + ensemble.fit$specificities
  
  error[i,1] <- expr.fit$results[1]
  error[i,2] <- meth.fit$results[1]
  error[i,3] <- cn.fit$results[1]
  error[i,4] <- ensemble.fit$results
  
  f1[i,1] <- expr.fit$f1score
  f1[i,2] <- meth.fit$f1score
  f1[i,3] <- cn.fit$f1score
  f1[i,4] <- ensemble.fit$f1score
  
  auc[i,1] <- expr.fit$auc
  auc[i,2] <- meth.fit$auc
  auc[i,3] <- cn.fit$auc
  auc[i,4] <- ensemble.fit$auc
  
  
  
}
#Aggregating avg information into one table
avg_error_expr <- mean(error[,1])
avg_error_meth <- mean(error[,2])
avg_error_cn <- mean(error[,3])
avg_error_ensemble <- mean(error[,4])
error_vec <- c(avg_error_expr,avg_error_meth,avg_error_cn,avg_error_ensemble)

avg_f1_expr <- mean(f1[,1])
avg_f1_meth <- mean(f1[,2])
avg_f1_cn <- mean(f1[,3])
avg_f1_ensemble <- mean(f1[,4])
f1_vec <- c(avg_f1_expr,avg_f1_meth,avg_f1_cn,avg_f1_ensemble)

avg_auc_expr <- mean(auc[,1])
avg_auc_meth <- mean(auc[,2])
avg_auc_cn <- mean(auc[,3])
avg_auc_ensemble <- mean(auc[,4])
auc_vec <- c(avg_auc_expr,avg_auc_meth,avg_auc_cn,avg_auc_ensemble)

information_matrix <- matrix(c(error_vec,f1_vec,auc_vec), nrow=3,ncol=4,byrow=T)
colnames(information_matrix) = c("Expr","Meth","CN", "Ensemble")
rownames(information_matrix) = c("Avg Error","Avg F1", "Avg AUC")
save(information_matrix, file=paste0("EnsembleOutput",1014,".RData"))

###PLOTTING THE ROC CURVES
avg_sens_expr <- expr_roc[,1]/5
avg_spec_expr <- expr_roc[,2]/5
df_expr <- data.frame(FPR = (1-avg_spec_expr), TPR = avg_sens_expr)

avg_sens_meth <- meth_roc[,1]/5
avg_spec_meth <- meth_roc[,2]/5
df_meth <- data.frame(FPR = (1-avg_spec_meth), TPR = avg_sens_meth)

avg_sens_cn <- cn_roc[,1]/5
avg_spec_cn <- cn_roc[,2]/5
df_cn <- data.frame(FPR = (1-avg_spec_cn), TPR = avg_sens_cn)

avg_sens_ensemble <- ensemble_roc[,1]/5
avg_spec_ensemble <- ensemble_roc[,2]/5
df_ensemble <- data.frame(FPR = (1-avg_spec_ensemble), TPR = avg_sens_ensemble)



rocs <- ggplot() +
  geom_line(data = df_expr, aes(x = FPR, y = TPR), color = "blue", size = 1.0, alpha = 0.7) + 
  geom_line(data = df_meth, aes(x = FPR, y = TPR), color = "green", size = 1.0, alpha = 0.7) + 
  geom_line(data = df_cn, aes(x = FPR, y = TPR), color = "red", size = 1.0, alpha = 0.7) +
  geom_line(data = df_ensemble, aes(x = FPR, y = TPR), color = "purple", size = 1.0, alpha = 0.7) +
  geom_abline(slope=1,intercept=0,linetype='dashed',alpha=0.7,color='grey')
  labs(title = "ROC Curves")

ggsave(rocs, paste0("rocs",1014,".png"))

 













