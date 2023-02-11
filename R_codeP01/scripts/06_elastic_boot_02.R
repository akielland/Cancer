#################################################################
##  Bootstrap  Elastic Net - using original sample as test set ##
#################################################################

## Elastic Net using bootstrap sample to train 1000 models
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## Test against original data sample
## output: proliferation.score correlation; SEM; Coefficient of genes

library(glmnet)


elastic_bootstrap_sample <- function(X, Y, ext.X, ext.Y, lambda.min=TRUE, method="pearson"){
  # one bootstrap sample
  # output: 
  # 1. correlation between prediction and full input sample  (pearson or spearman should be set)
  # 2. Coefficients of variables 
  # 3. MSE
  n = length(Y)
  
  int <- sample.int(n, size = n, replace = TRUE)
  X_train = X[int,]
  Y_train = Y[int]
  
  alphas = seq(0, 1, by=0.1)
  
  # To store results
  fit.cv <- list()
  cvm = rep(0, length(alphas))
  
  # run elastic net with cross-validation and alpha search
  for (i in 1:length(alphas)){
    fit.cv[[i]] <- cv.glmnet(X_train, Y_train, alpha = alphas[i], nfolds = 5)
    cvm[i] <- min(fit.cv[[i]]$cvm)
  }
  
  opt.idx <- which.min(cvm)
  fit.cv <- fit.cv[[opt.idx]]
  
  co <- as.vector(coef(fit.cv, s="lambda.min"))
  
  if (lambda.min == TRUE) {
    pred_values = predict(fit.cv, newx = ext.X, type = "response", s = "lambda.min")      
  }
  if (lambda.min == FALSE) {
    pred_values = predict(fit.cv, newx = ext.X, type = "response", s = "lambda.1se")      
  }
  cor <- suppressWarnings(cor(pred_values, ext.Y, method = method))
  MSE <- mean((ext.Y - pred_values)^2)
  return(list(cor=cor, co=co, MSE=MSE))
}

elastic_bootstrap_sample(X, Y, X, Y, TRUE)$cor


# 1000 bootstrap fits
elastic_boot = function(X, Y, ext.X, ext.Y, n_bootstraps=1000, method = "pearson") {
  # run many bootstraps
  # output: - vector with correlations 
  #         - matrix with features coefficients
  #         - vector with SEM
  cor_vec <- rep(NA, n_bootstraps)
  co_matrix <- matrix(NA, nrow = dim(X)[2], ncol = n_bootstraps)
  rownames(co_matrix) <- colnames(X) 
  MSE_vec <- rep(NA, n_bootstraps)
  
  for (i in c(1:n_bootstraps)) {
    out <- elastic_bootstrap_sample(X, Y, ext.X, ext.Y, TRUE, method = method)
    # browser()
    cor_vec[i] <- out$cor
    co_matrix[,i] <- out$co[-1]
    MSE_vec[i] <- out$MSE
    cat(i, "")
  } 
  return(list(cor_vec=cor_vec, co_matrix=co_matrix, MSE_vec = MSE_vec))
}

# 1000 bootstrap fits
# RUN: e_b_obj_771_p
set.seed(123)
e_b_obj_771_p <- elastic_boot(X, Y, X, Y)
head(e_b_obj_771_p$co_matrix)[,1:6]
save(e_b_obj_771_p, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/e_b_obj_771_p.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/e_b_obj_771_p.RData")

# RUN: e_b_obj_771_ROR_p
set.seed(123)
e_b_obj_771_ROR_p <- elastic_boot(X, Y, X, Y)
head(e_b_obj_771_ROR_p$co_matrix)[,1:6]
save(e_b_obj_771_ROR_p, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/e_b_obj_771_ROR_p.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/e_b_obj_771_ROR_p.RData")









eb_object <- elastic_boot(X, Y, X, Y, 1000)

eb_object
histogram(eb_object$cor_vec, breaks = 99,
          xlab = "Proliferation score", 
          main = "Elastic net (trail 1, arm Letro+Ribo)")


# Various objects
save(eb_object, file="eb_object_6genes.RData")
save(eb_object, file="eb_object_Nodes.RData")
save(eb_object, file="eb_object_AllGenes_Prolif.RData")
save(eb_object, file="eb_object_ALLGenes_ROR.RData")

# stored object can be loaded to save time
load("eb_object_6Genes.RData")
load("eb_object_Nodes.RData")
load("eb_object_AllGenes_Prolif.RData")
load("eb_object_ALLGenes_ROR.RData")

## Summery result of lasso bootstrap object
# plot cross-validated error as a function of lambda and alpha

# Correlation
cor_vec <- as.numeric(eb_object$cor_vec)
sum(is.na(cor_vec))/length(cor_vec)    # fraction of NA
mean(cor_vec, na.rm=TRUE)
var(cor_vec, na.rm=TRUE)
par(mfrow=c(1,1))
hist(cor_vec, breaks = 100)

# MSE
MSE_vec <- eb_object$MSE_vec
mean(MSE_vec)
var(MSE_vec)
hist(MSE_vec, breaks=100)

## Features
genes_mean <- rowMeans(eb_object$co_matrix)
apply(eb_object$co_matrix, 1, var)
sort(genes_mean, decreasing = TRUE)
sort(abs(genes_mean), decreasing = TRUE)[1:3]

# Bar diagram on features ordered by amount of times selected
df_ <- as.data.frame(genes_mean)
df_["genes_names"] <- row.names(df_)
df_ <- df_[order(df_$genes_mean, decreasing = TRUE),]
df_pluss <- df_[1:2,]
gene_n <- dim(df_)[1];   gene_n2 <- gene_n-1
df_minus <- df_[gene_n2: gene_n, ]

# lock in factor level order
df_$genes_names <- factor(df_$genes_names, levels = df_$genes_names)

# plot
ggplot(data=df_, aes(x=genes_names, y=genes_mean)) + 
  geom_bar(stat="identity") + coord_flip()

ggplot(data=df_, aes(x=genes_names, y=genes_mean)) + 
  geom_col() + coord_flip()




