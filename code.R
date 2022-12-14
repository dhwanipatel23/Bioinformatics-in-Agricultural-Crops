# setting the seed value to get consistent results
set.seed(50)

# load the necessary libraries

library(GEOquery)
library(limma)
library(umap)
library(Biobase)
library(dplyr)
library(caret)
library(randomForest)
library(genefilter)
library(tuneRanger)
library(mlr)
library(e1071)
library(tidyverse) 
library(rpart)
library(rpart.plot)
library(pROC)


# loading the dataset

gset <- getGEO("GSE50770", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8672", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

## Data preprocessing

# Log transformation of the dataset
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

## Exploratory data analysis

# Visualization of data using the R-script of GEO2R  

fvarLabels(gset) <- make.names(fvarLabels(gset)) # make proper column names to match top table 


gsms <- "00022211133344455" # group membership for all samples
sml <- strsplit(gsms, split="")[[1]]

# assign samples to groups and set up design matrix

gs <- factor(sml)
groups <- make.names(c("normal","cold","acidic","drought","saline","alkaline"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients

cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names

ct <- 2        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

ct <- 4        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# General expression data analysis
exa <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE50770", "/", annotation(gset), sep ="")
boxplot(exa[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE50770", "/", annotation(gset), " value distribution", sep ="")
plotDensities(exa, group=gs, main=title, legend ="topright")


# Feature selection:

# Filtering genes to reduce dimension
ffun <- filterfun(pOverA(0.20,100), cv(0.9,5))
t.fil <- genefilter(ex,ffun) 

# apply filter, and put expression back on log scale
small.eset <- log2(ex[t.fil,])


dim(ex) # 24132    17

dim(small.eset) # 4610   17

# Transposing the data to get genes as features
filtered.data <- t(small.eset)

dim(filtered.data) # 17 4610


# defining class vector
labels <- c('1','1','1','2','2','2','3','3','3','4','4','4','5','5','5','6','6')

###################

# performing leave-one-out cross-validation on the data

# creating empty vectors to store predicted labels and actual labels inside the loop
rfc.label <- c()
svmr.label <- c()
svml.label <- c()
dt.label <- c()
actual.label <- c()

# defining index for loocv loop
data.idx <- c(1:nrow(filtered.data))
for (ts in c(1:length(data.idx))){
  set.seed(50)
  
  # train index and test index
  train.idx <- setdiff(data.idx, ts)
  test.idx <- ts
  
  # generating train data and labels as per the train index
  data.tr <- filtered.data[train.idx,]
  label.tr <- labels[train.idx]
  
  print('Train labels: ',label.tr)
  
  # generating test data and labels as per the test index
  data.ts <- t(filtered.data[test.idx,])
  label.ts <- labels[test.idx]
  print('Test labels: ',label.ts)
  
  tr.df<-as.data.frame(data.tr)
  ts.df<-as.data.frame(data.ts)

  
  # Hyper-parameter tuning for Random Forest
  # storing train data to an empty dataframe  
  df <- data.frame()
  df<-data.frame(data.tr)
  df$dflabel <- label.tr
  
  
  rf.task = makeClassifTask(data = df, target = "dflabel")
  
  # getting the best parameters for the model using tuneRanger
  
  res = tuneRanger(rf.task, measure = list(multiclass.brier), num.trees = 700, 
                   num.threads = 2, iters = 20)
  
  param <- res$recommended.pars
  
  # Random Forest Classifier with tuned parameters
  rf_model = randomForest(x = data.tr, 
                     y = factor(label.tr), 
                     num.threads = 2,
                     verbose = FALSE,
                     respect.unordered.factors = order, 
                     num.trees = 700, 
                     mtry = param$mtry,
                     min.node.size = param$min.node.size,
                     sample.fraction = param$sample.fraction,
                     replace = FALSE)
  
  # SVM model with radial kernel
  svm_model.r = svm(label.tr ~., data=data.tr, 
                   type ="C-classification", kernal="radial", 
                   gamma=0.1, cost=10)
  
  # SVM model with linear kernel
  svm_model.l = svm(label.tr ~., data=data.tr, 
                   type="C-classification", kernal="linear", 
                   gamma=0.1, cost=10)
  
  # Decision tree model
  dt_model <- rpart(label.tr ~., data = tr.df, 
                    method = 'class',minsplit = 10, minbucket=3)
  
  

  # predicting the class with each model
  rfc_pred = predict(rf_model, newdata = data.ts, type = 'response')
  svmr_pred = predict(svm_model.r, newdata = data.ts, type = 'response')
  svml_pred = predict(svm_model.l, newdata = data.ts, type = 'response')
  dt_pred = predict(dt_model, newdata = ts.df, type = 'class')
  
  
  
  # updating the predicted label and actual label in the empty vectors
  rfc.label <- append(rfc.label,rfc_pred)
  svmr.label <- append(svmr.label,svmr_pred)
  svml.label <- append(svml.label,svml_pred)
  dt.label <- append(dt.label,dt_pred)
  
  actual.label <- append(actual.label,label.ts)
  
  
}

# Display matrix of labels and determine accuracy using diagonal values
rfc.table <- table(rfc.label,actual.label)
svmr.table <- table(svmr.label,actual.label)
svml.table <- table(svml.label,actual.label)
dt.table <- table(dt.label,actual.label)

rfc.accuracy <- sum(diag(rfc.table))/sum(rfc.table)
svmr.accuracy <- sum(diag(svmr.table))/sum(svmr.table)
svml.accuracy <- sum(diag(svml.table))/sum(svml.table)
dt.accuracy <- sum(diag(dt.table))/sum(dt.table)

print(rfc.table)
print(rfc.accuracy) 

print(svmr.table)
print(svmr.accuracy) 

print(svml.table)
print(svml.accuracy) 

print(dt.table)
print(dt.accuracy) 


# Get roc curve for random forest

result <- pROC::multiclass.roc(as.numeric(rfc.label), as.numeric(actual.label))

plot.roc(result$rocs[[1]], print.auc=T, legacy.axes = T, print.auc.adj = c(0,1))

plot.roc(result$rocs[[2]], add=T, col = 'red', print.auc = T, legacy.axes = T, print.auc.adj = c(0,3))

plot.roc(result$rocs[[3]],add=T, col = 'blue', print.auc=T, legacy.axes = T, print.auc.adj = c(0,5))

plot.roc(result$rocs[[4]],add=T, col = 'orange', print.auc=T, legacy.axes = T, print.auc.adj = c(0,7))

plot.roc(result$rocs[[5]],add=T, col = 'brown', print.auc=T, legacy.axes = T, print.auc.adj = c(0,9))

plot.roc(result$rocs[[6]],add=T, col = 'pink', print.auc=T, legacy.axes = T, print.auc.adj = c(0,11))

