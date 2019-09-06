# packages:
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

# loading data
screening <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/screening.rds")
meth <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/methylation.rds")
copynum <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/copynumber.rds")
expression <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/expression.rds")

# copy number data based on significant genes for each drug (15 drugs btwn 20% and 80% EFFECT-INEFFECTIVE ratio)
copynum.1016 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1016_0.1.rds")
copynum.1007 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1007_0.1.rds")
copynum.1026 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1026_0.1.rds")
copynum.1006 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1006_0.1.rds")
copynum.1058 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1058_0.1.rds")
copynum.1015 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1015_0.1.rds")
copynum.1014 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1014_0.1.rds")
copynum.1001 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1001_0.1.rds")
copynum.1053 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1053_0.1.rds")
copynum.1060 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1060_0.1.rds")
copynum.1011 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1011_0.1.rds")
copynum.1054 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1054_0.1.rds")
copynum.1037 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1037_0.1.rds")
copynum.1066 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1066_0.1.rds")
copynum.1008 <- readRDS("C:/Users/yasha/Desktop/BDSI/Data Mining/Project/Copy Number/copynum1008_0.1.rds")

# data imputation: use copynum.1016
View(colSums(copynum.1016)) # NAs exist; thus, we must impute

## MICE package: Assumptions are that probabilty that a value is missing depends 
# on only observed values and thus can be predicted using sed variables. Imputation 
# is done on a variable by variable basis specifiying imputation model used.
install.packages("mice")
library(mice)
View(md.pattern(copynum.1016))

temp.copynum.1016 <- mice(copynum.1016, method = "pmm", seed = 500)








