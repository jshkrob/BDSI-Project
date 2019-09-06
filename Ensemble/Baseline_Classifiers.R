
# Packages:
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

cell_line_type <- read_excel("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Cell_Lines_Details.xlsx", sheet = "Cell line details" )
colnames(cell_line_type) <- c("CELL_LINE_NAME", "CL", "WES", "CNA", "GENE EXP", "METHY",
                              "Drug Response", "Site", "Tissue Site 2", "Cancer Type", 
                              "MSI", "Screen Medium", "Growth Properties")


cell_line_type <- cell_line_type[, c("CELL_LINE_NAME", "Site")]
screening <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/screening.rds")
expression <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/expression.rds")
meth <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/methylation.rds")
copynum <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/copynumber.rds")

screening <- screening %>% left_join(cell_line_type, by = "CELL_LINE_NAME")
  
# drug-specific data:
effective.1054 <- subset(screening, DRUG_ID_lib == "1016")[,c("CL","EFFECT", "Site")]
expr.1054 <- expression[as.character(effective.1054$CL), ]

meth.1054 <- meth[as.character(effective.1054$CL), ]
meth.1054 <- meth.1054[ , !is.na(colSums(meth))]
copynum.1054 <- copynum[as.character(effective.1054$CL), ]

## First, create a one-hot encoding for given drug: colnames = site.{tissue}
indicator_site <- data.frame(to.dummy(effective.1054$Site, prefix = "site"), 
                             row.names = effective.1054$CL)
## perform classifiers:
ERROR_site <- data.frame(nb = numeric(10), glm = numeric(10), lasso = numeric(10), ridge = numeric(10), elas = numeric(10))


## Now subset (for the given drug) each data set (meth, expr, copynum) by the type of tissue, and create a "centered"
## data frame for each data set based on the tissue feature means:

expr.1054$Site <- effective.1054$Site # append Site column to group by tissues...
meth.1054$Site <- effective.1054$Site # append Site column to group by tissues...
copynum.1054$Site <- effective.1054$Site # append Site column to group by tissues...

scale_this <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  return(x - mu)
}

############################################ Centering of data:

# Indicator Matrix:
indicator_site <- indicator_site %>% 
  add_column(Site = effective.1054$Site) %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  ungroup() %>% 
  dplyr::select(-Site)

# EXPR
expr.1054_c <- expr.1054 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)

#rownames(expr.1054_c) <- rownames(expr.1054)

# METH
meth.1054_c <- meth.1054 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)

#rownames(meth.1054_c) <- rownames(expr.1054)

# COPY
copynum.1054_c <- copynum.1054 %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  transmute_at(vars(-Site), scale_this) %>% 
  ungroup() %>% 
  dplyr::select(-Site)

#rownames(copynum.1054_c) <- rownames(expr.1054)

# accordingly, do the same with the effect (to remove infreqeunt tissue cases...)

Y <- effective.1054[, c("EFFECT", "Site")]
Y_c <- Y %>% 
  group_by(Site) %>% 
  filter(n() >= 5) %>% 
  ungroup() %>% 
  dplyr::select(-Site) %>% 
  pull(EFFECT)

## output: {dataset}.{drugid}_c and Y_c 


# choose best classification method to use using the indicator matrix "indicator_c" and effect classification "Y_c" 

for(i in 1:10) { 
  # choosing training and test
  n_obs <- nrow(indicator_site)
  train.ix = sample(1:n_obs, round(n_obs* 0.6))
  ytrain <- Y_c[train.ix] 
  xtrain <- indicator_site[train.ix,]
  ytest <- Y_c[-train.ix]
  xtest <- indicator_site[-train.ix, ]
  
  # naive bayes:
  model_nb = naive_bayes(x = xtrain, y = as.factor(ytrain))
  ytrainhat = predict(model_nb, xtrain)
  ytesthat = predict(model_nb, xtest)
  #train_ERROR_histog = 1 - mean(ytrain == ytrainhat)
  testerror = 1 - mean(ytest == ytesthat)
  ERROR_site[i,1] <- testerror  
  
  # Logistic Regression:
  train <- data.frame(ytrain, xtrain)
  model_glm <- glm(ytrain ~ ., family = binomial(link = "logit"), data = train) # increasing allowed iterations beyond default of 25...
  ytrainhat <- ifelse(model_glm$fitted.values >= 0.5, 1, 0)
  ytesthat <- ifelse(predict(model_glm, xtest, type = "response") >= 0.5, TRUE, FALSE)
  #trainingerror = 1 - mean(ytrain==ytrainhat)
  testerror = 1 - mean(ytest==ytesthat)
  ERROR_site[i,2] <- testerror
  
  fit <- cv.glmnet(as.matrix(indicator_site), y = Y_c, family = "binomial", alpha = 1)
  lambda_ind <- fit$lasso.min
  
  # Lasso
  model_lasso <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", alpha = 1, lambda = lambda_ind)
  ytrainhat <- ifelse(predict(model_lasso, as.matrix(xtrain), type = "response") >= 0.5, TRUE, FALSE)
  ytesthat <- ifelse(predict(model_lasso, as.matrix(xtest), type = "response") >= 0.5, TRUE, FALSE)
  #trainingerror <- mean(ytrainhat != ytrain)
  testerror <- mean(ytesthat != ytest)
  ERROR_site[i,3] <- testerror  
  
  fit <- cv.glmnet(as.matrix(indicator_site), y = Y_c, family = "binomial", alpha = 0)
  lambda_ind <- fit$lasso.min
  
  # Ridge
  model_ridge <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", alpha = 0, lambda = lambda_ind)
  ytrainhat <- ifelse(predict(model_ridge, as.matrix(xtrain), type = "response") >= 0.5, TRUE, FALSE)
  ytesthat <- ifelse(predict(model_ridge, as.matrix(xtest), type = "response") >= 0.5, TRUE, FALSE)
  #trainingerror <- mean(ytrainhat != ytrain)
  testerror <- mean(ytesthat != ytest)
  ERROR_site[i,4] <- testerror  
  
  
  fit <- cv.glmnet(as.matrix(indicator_site), y = Y_c, family = "binomial", alpha = 0.5)
  lambda_ind <- fit$lasso.min
  
  # Elastic Net:
  model_elas <- glmnet(x = as.matrix(xtrain), y = ytrain, family = "binomial", alpha = 0.5, lambda = lambda_ind)
  ytrainhat <- ifelse(predict(model_elas, as.matrix(xtrain), type = "response") >= 0.5, TRUE, FALSE)
  ytesthat <- ifelse(predict(model_elas, as.matrix(xtest), type = "response") >= 0.5, TRUE, FALSE)
  #trainingerror <- mean(ytrainhat != ytrain)
  testerror <- mean(ytesthat != ytest)
  ERROR_site[i,5] <- testerror
}  

# glm for "site" and "histography"
summary(ERROR_site)
