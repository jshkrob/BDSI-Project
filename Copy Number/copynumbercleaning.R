# Final Project
screening = readRDS("C:/Users/Corri Sept/Documents/R/LearningR/screening.rds")
methylation = readRDS("C:/Users/Corri Sept/Documents/R/LearningR/methylation.rds")
expr = readRDS("C:/Users/Corri Sept/Documents/R/LearningR/expression.rds")
copynum = readRDS("C:/Users/Corri Sept/Documents/R/LearningR/copynumber.rds")
effective <- subset(screening)[,c("CL","CELL_LINE_NAME","EFFECT")]

dim(copynum)
View(screening)
head(copynum)
View(copynum)

# Dealing with copy number variation data...so how many copies of a gene per patient
# I would like to do feature selection followed by clustering

# There is missing data...need to figure out if those variables are important or not (in terms of drug effectiveness)
# Specifically, want to determine which genes are important for predicting whether or not a drug 
# is effective for a particular cell line. 
# Therefore, need to check everything 59 times. 

# First, create a vector of drugs that I want to loop
# through, excluding ones that are completely ineffective or completely effective

library(tidyverse)

# Look at what percentage of the time drugs are effective
View(screening %>% group_by(DRUG_ID_lib) %>% summarize(proportion = mean(EFFECT))
     %>% arrange(desc(proportion)))

# Make variable of proportion of proportion of time drugs are EFFECTIVE
perceffect <- screening %>% group_by(DRUG_ID_lib) %>% summarize(proportion = mean(EFFECT)) %>% 
      arrange(desc(proportion))

View(perceffect)

# Get drug IDs of ones that are effective between 20% and 80% of the time
indices <- which(perceffect$proportion<=0.8 & perceffect$proportion>=0.2,arr.ind=T)
drugsclassify <- perceffect[indices, 1] # drugsclassify is a list of all drug IDs we are using

# Just wanted to check how many observations for this one, since Ritoban mentioned it had less than 150
# It has only 115, but we discussed (7/3/2019) and have decided to retain it
effective.1037 <- subset(screening, DRUG_ID_lib == "1037")[,c("CL","CELL_LINE_NAME","EFFECT")]

# Note that drugclassify includes 15 drug ID names...this is significantly lower than our beginning 59
# Also note that 8 drugs are effective 0% of the time...concerning
# Now subset the other datasets using drugclassify, which is a vector of all drugs we are using

screening15drugs <- subset(screening, DRUG_ID_lib %in% drugsclassify$DRUG_ID_lib)[,c("DRUG_ID_lib","CL","CELL_LINE_NAME","EFFECT")]

# Make sure subsetted correctly
checkscreening <- unique(screening15drugs$DRUG_ID_lib)
View(checkscreening) # looks good, we've got 15 drugs and they're the right ones

# Now I'm going to save my screening data for the specific drugs
saveRDS(screening15drugs,"C:/Users/Corri Sept/Documents/R/LearningR/screening15drugs.rds")


# Function for t-test
get.p <- function(dat, labs){
  # split the data into effective and ineffective
  effect <- dat[labs]
  ineffect <- dat[!labs]
  
  # # calculate the two sample means
  # effect.bar <- mean(effect, na.rm=TRUE)
  # ineffect.bar <- mean(ineffect, na.rm=TRUE)
  # 
  # # calculate the two sample variances
  # v.effect <- var(effect, na.rm=TRUE)
  # v.ineffect <- var(ineffect, na.rm=TRUE)
  # 
  # # calculate the sample sizes
  # n.effect <- length(effect)
  # n.ineffect <- length(ineffect)
  # 
  # # calculate the sd
  # s_pooled = (n.effect*v.effect +n.ineffect*v.ineffect)/(n.effect+n.ineffect-2)
  # s <- sqrt(s_pooled*(1/n.effect + 1/n.ineffect))
  # #s <- sqrt((v.effect/n.effect) + (v.ineffect/n.ineffect))
  # 
  # # calculate the test statistic
  # T_abs <- abs((effect.bar - ineffect.bar)/s)
  # pval = 2*(1-pt(T_abs,df = n.effect+n.ineffect - 2))
  
  pval = t.test(effect,ineffect,alternative = "two.sided", na.action = na.omit)$p.value
  return(pval)
}

# Specifically, want to determine which genes are important for predicting whether or not a drug 
# is effective for a particular cell line. 
# Therefore, need to check everything 59 times. First, create a vector of drugs that I want to loop
# through, excluding ones that are completely ineffective or completely effective

# Remove columns w missing values for every observation
# remove columns with missing values and whose t-test aren't significant
# impute remaining ones with NAs

# FIRST, remove columns with missing values for every observation
# Percentage of missing values (across cell lines) for each gene
df$percentNA <- data.frame(colMeans(is.na(copynum)))

percentNA <- colMeans(is.na(copynum))
View(percentNA)

# Extract row names of genes where all values are missing
#allna <- which(df$percentNA==1,arr.ind=T)
#rowname <- row.names(df[allna[,1],] ) #row names

allna <- which(percentNA==1,arr.ind=T)
rownamesubset <- percentNA[allna]
rowname <- names(rownamesubset)

library(dplyr)

# Genes that don't have all missing values
somemiss <- setdiff(colnames(copynum), rowname)

# Remove (gene) columns with missing values for every observation
copynumsomemiss <- copynum[,somemiss]
dim(copynumsomemiss) # It's the right dimension, looks like it worked!

View(copynumsomemiss) # Looks like it kept the correct columns in
checkNA <- colMeans(is.na(copynumsomemiss)) # No longer any genes with all missing, and dimensions check out, so looks like we're good
View(checkNA)

# Now save the dataset with genes that DO NOT have ONLY missing values
saveRDS(copynumsomemiss,"C:/Users/Corri Sept/Documents/R/LearningR/copynumsomemiss.rds")

# First, I'm going to apply this function for one treatment. I'll use 1014 
effective.1014 <- subset(screening, DRUG_ID_lib == "1014")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1014 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1014$CL)

# Create empty dataframe for p-values
pvals <- data.frame(cell_lines = colnames(copynumsomemiss))

pvals$p_1014 = apply(copynum.1014,2,get.p, effective.1014$EFFECT)

View(pvals %>% arrange(desc(p)))



# Now I'm going to apply this for all 15 treatments

# Make empty data frame of p-values
pvals <- data.frame(genes = colnames(copynumsomemiss))
View(pvals)


# I'm sure there's a more elegant way to do this, but I'll just brute force it and subset / apply function for all 15 drugs

### 1007
# Subset 
effective.1007 <- subset(screening15drugs, DRUG_ID_lib == "1007")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1007 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1007$CL) # making sure have same cell lines
# Apply function
pvals$p_1007 = apply(copynum.1007,2,get.p, effective.1007$EFFECT)

### 1026
# Subset
effective.1026 <- subset(screening15drugs, DRUG_ID_lib == "1026")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1026 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1026$CL)
# Apply function
pvals$p_1026 = apply(copynum.1026,2,get.p, effective.1026$EFFECT)

### 1016
# Subset
effective.1016 <- subset(screening15drugs, DRUG_ID_lib == "1016")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1016 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1016$CL)
# Apply function
pvals$p_1016 = apply(copynum.1016,2,get.p, effective.1016$EFFECT)
# GOT AN ERROR... NOT ENOUGH Y OBSERVATIONS. MANUALLY CHECK

effect1016sub <- copynum.1016[effective.1016$EFFECT]
ineffect1016sub <- copynum.1016[!effective.1016$EFFECT]

# This is just to check to make sure that loop works correctly
for(i in 1:5) {
  get.p(copynum.1016[,i],effective.1016$EFFECT)
  i = i+1
}

# Figure out which gene is giving me the error (note that below code was after finding
# 5823 and 5824; I wanted to check genes that hadn't been checked yet)
for(i in 5825:46157) {
  get.p(copynum.1016[,i],effective.1016$EFFECT)
  i = i+1
}
# stop at i = 5823 
# stop at i = 5824
get.p(copynum.1016[,5823],effective.1016$EFFECT)

get.p(copynum.1016[,5824],effective.1016$EFFECT)

# As it turns out, 5823 and 5824 are giving me errors. I need to remove these columns

# I'm going to do this cautiously because I don't want to screw up and remove everything
dim(copynum.1016) # 233 x 46157
copynum1016rc <- copynum.1016[,-5823]

dim(copynum1016rc) # Correct dimension, 233 x 46156
copynum1016rc2 <- copynum1016rc[,-5823] # because I just removed column 5823, column 5824 has become 5823
dim(copynum1016rc2) # Correct dimension, 233 x 46155


# Make sure this worked
for(i in 1:46155) {
  get.p(copynum1016rc2[,i],effective.1016$EFFECT)
  i = i+1
}


# I'll run the t-test now. I'll know for sure by doing this if I've correctly removed columns / 
# If I've removed the right ones
# NOTE that I am not removing columns from effective.1016...this is because the problem is I 
# have missing values in copy numbers for genes. effective.1016 does not have this information, 
# only copynum.1016 does.

# Note that r studio is throwing me errors because my p-value dataframe is 46157 rows and because
# I deleted troublesome genes, I only have 46155 for this particular drug. It gets its own dataframe
pvals_1016 <- data.frame(genes = colnames(copynum1016rc2))

# Apply function
pvals_1016$p_1016 = apply(copynum1016rc2,2,get.p, effective.1016$EFFECT)

# Alright, seems to be working now. Checked dimension and viewed the dataset itself.

### 1006
# Subset
effective.1006 <- subset(screening15drugs, DRUG_ID_lib == "1006")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1006 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1006$CL)
# Apply function
pvals$p_1006 = apply(copynum.1006,2,get.p, effective.1006$EFFECT)

### 1058
# Subset
effective.1058 <- subset(screening15drugs, DRUG_ID_lib == "1058")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1058 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1058$CL)
# Apply function
pvals$p_1058 = apply(copynum.1058,2,get.p, effective.1058$EFFECT)

### 1015
# Subset
effective.1015 <- subset(screening15drugs, DRUG_ID_lib == "1015")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1015 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1015$CL)
# Apply function
pvals$p_1015 = apply(copynum.1015,2,get.p, effective.1015$EFFECT)

### 1014
# Subset
effective.1014 <- subset(screening15drugs, DRUG_ID_lib == "1014")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1014 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1014$CL)
# Apply function
pvals$p_1014 = apply(copynum.1014,2,get.p, effective.1014$EFFECT)

### 1001
# Subset
effective.1001 <- subset(screening15drugs, DRUG_ID_lib == "1001")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1001 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1001$CL)
# Apply function
pvals$p_1001 = apply(copynum.1001,2,get.p, effective.1001$EFFECT)

### 1053
# Subset
effective.1053 <- subset(screening15drugs, DRUG_ID_lib == "1053")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1053 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1053$CL)
# Apply function
pvals$p_1053 = apply(copynum.1053,2,get.p, effective.1053$EFFECT)

### 1060
# Subset
effective.1060 <- subset(screening15drugs, DRUG_ID_lib == "1060")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1060 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1060$CL)
# Apply function
pvals$p_1060 = apply(copynum.1060,2,get.p, effective.1060$EFFECT)

### 1011
# Subset
effective.1011 <- subset(screening15drugs, DRUG_ID_lib == "1011")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1011 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1011$CL)
# Apply function
pvals$p_1011 = apply(copynum.1011,2,get.p, effective.1011$EFFECT)

### 1054
# Subset
effective.1054 <- subset(screening15drugs, DRUG_ID_lib == "1054")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1054 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1054$CL)
# Apply function
pvals$p_1054 = apply(copynum.1054,2,get.p, effective.1054$EFFECT)

### 1037
# Subset
effective.1037 <- subset(screening15drugs, DRUG_ID_lib == "1037")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1037 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1037$CL)
# Apply function
pvals$p_1037 = apply(copynum.1037,2,get.p, effective.1037$EFFECT)

### 1066
# Subset
effective.1066 <- subset(screening15drugs, DRUG_ID_lib == "1066")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1066 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1066$CL)
# Apply function
pvals$p_1066 = apply(copynum.1066,2,get.p, effective.1066$EFFECT)

### 1008
# Subset
effective.1008 <- subset(screening15drugs, DRUG_ID_lib == "1008")[,c("CL","CELL_LINE_NAME","EFFECT")]
copynum.1008 <- subset(copynumsomemiss, rownames(copynumsomemiss) %in% effective.1008$CL)
# Apply function
pvals$p_1008 = apply(copynum.1008,2,get.p, effective.1008$EFFECT)


# p-values from t-tests for every drug except for 1016
saveRDS(pvals,"C:/Users/Corri Sept/Documents/R/LearningR/copynum_p_values.rds")

# p-value for drug 1016
saveRDS(pvals_1016,"C:/Users/Corri Sept/Documents/R/LearningR/copynum_p_values_1016.rds")

# NEXT, I'm going to subset so that anything with a value less than .15 (FOR THAT PARTICULAR TREATMENT) is cut out. 
# By subset, I'm including genes that don't have any NAs. I'm just doing general feature selection, and I'll have 15 datasets.
# Following, I'm going to check to see if I still have genes with NAs. If I do, then I'll impute those genes for those treatments.

# Number of genes kept in is commented to the right of the subsetting
# ALSO note that I'm subsetting from datasets that were (a) made from copynumsomemiss, the one with no genes that are all missing, and 
# (b) that have cell lines that match up with the tested ones (so ones that were actually in the asssay)

# 1016
p1016sig <- subset(pvals_1016, p_1016<0.1) # 24317
copynum1016_0.1 <- copynum1016rc2[,colnames(copynum1016rc2) %in% p1016sig$genes] # 233 x 24317
saveRDS(copynum1016_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1016_0.1.rds")

# 1007
p1007sig <- subset(pvals, p_1007<0.1, select = c(genes, p_1007)) # 14794
copynum1007_0.1 <- copynum.1007[,colnames(copynum.1007) %in% p1007sig$genes] # 246 x 14794
saveRDS(copynum1007_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1007_0.1.rds")

# 1026
p1026sig <- subset(pvals, p_1026<0.1, select = c(genes, p_1026)) # 6721
copynum1026_0.1 <- copynum.1026[,colnames(copynum.1026) %in% p1026sig$genes] # 306 x 6721
saveRDS(copynum1026_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1026_0.1.rds")

# 1006
p1006sig <- subset(pvals, p_1006<0.1, select = c(genes, p_1006)) # 22316
copynum1006_0.1 <- copynum.1006[,colnames(copynum.1006) %in% p1006sig$genes] # 269 x 22316
saveRDS(copynum1006_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1006_0.1.rds")

# 1058
p1058sig <- subset(pvals, p_1058<0.1, select = c(genes, p_1058)) # 4801
copynum1058_0.1 <- copynum.1058[,colnames(copynum.1058) %in% p1058sig$genes] # 191 x 4801
saveRDS(copynum1058_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1058_0.1.rds")

# 1015
p1015sig <- subset(pvals, p_1015<0.1, select = c(genes, p_1015)) # 22,209
copynum1015_0.1 <- copynum.1015[,colnames(copynum.1015) %in% p1015sig$genes] # 166 x 22,209
saveRDS(copynum1015_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1015_0.1.rds")

# 1014
p1014sig <- subset(pvals, p_1014<0.1, select = c(genes, p_1014)) # 11,994
copynum1014_0.1 <- copynum.1014[,colnames(copynum.1014) %in% p1014sig$genes] # 246 x 11,994
saveRDS(copynum1014_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1014_0.1.rds")

# 1001
p1001sig <- subset(pvals, p_1001<0.1, select = c(genes, p_1001)) # 36,517
copynum1001_0.1 <- copynum.1001[,colnames(copynum.1001) %in% p1001sig$genes] # 152 x 36,517
saveRDS(copynum1001_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1001_0.1.rds")

# 1053
p1053sig <- subset(pvals, p_1053<0.1, select = c(genes, p_1053)) # 4317
copynum1053_0.1 <- copynum.1053[,colnames(copynum.1053) %in% p1053sig$genes] # 170 x 4317
saveRDS(copynum1053_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1053_0.1.rds")

# 1060
p1060sig <- subset(pvals, p_1060<0.1, select = c(genes, p_1060)) # 7498
copynum1060_0.1 <- copynum.1060[,colnames(copynum.1060) %in% p1060sig$genes] # 251 x 7498
saveRDS(copynum1060_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1060_0.1.rds")

# 1011
p1011sig <- subset(pvals, p_1011<0.1, select = c(genes, p_1011)) # 39819
copynum1011_0.1 <- copynum.1011[,colnames(copynum.1011) %in% p1011sig$genes] # 328 x 39,819
saveRDS(copynum1011_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1011_0.1.rds")

# 1054
p1054sig <- subset(pvals, p_1054<0.1, select = c(genes, p_1054)) # 22061
copynum1054_0.1 <- copynum.1054[,colnames(copynum.1054) %in% p1054sig$genes] # 243 x 22,061
saveRDS(copynum1054_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1054_0.1.rds")

# 1037
p1037sig <- subset(pvals, p_1037<0.1, select = c(genes, p_1037)) # 21624
copynum1037_0.1 <- copynum.1037[,colnames(copynum.1037) %in% p1037sig$genes] # 115 x 21624
saveRDS(copynum1037_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1037_0.1.rds")

# 1066
p1066sig <- subset(pvals, p_1066<0.1, select = c(genes, p_1066)) # 2897
copynum1066_0.1 <- copynum.1066[,colnames(copynum.1066) %in% p1066sig$genes] # 179 x 2897
saveRDS(copynum1066_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1066_0.1.rds")

# 1008
p1008sig <- subset(pvals, p_1008<0.1, select = c(genes, p_1008)) # 28,852
copynum1008_0.1 <- copynum.1008[,colnames(copynum.1008) %in% p1008sig$genes] # 364 x 28,852
saveRDS(copynum1008_0.1,"C:/Users/Corri Sept/Documents/R/LearningR/copynum1008_0.1.rds")
