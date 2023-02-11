## Ridge using bootstrap sample to train 1000 models
## Cross-validation (5-fold) used for training/tuning lambda and selecting features
## Test against original data sample and leave-out bootstrap sample
## output: proliferation.score correlation; MSE; coefficients of best model in each bootstrap

library(glmnet)

##########################################################
##  Bootstrap  Ridge - using original sample as test set ##
##########################################################


ridge_bootstrap_sample <- function(X, Y, ext.X, ext.Y, lambda.min=TRUE, method="pearson"){
  # one bootstrap sample
  # output: 
  # 1. correlation between prediction and full input sample  (pearson or spearman should be set)
  # 2. Coefficients of variables 
  # 3. MSE
  n = length(Y)
  
  int <- sample.int(n, size = n, replace = TRUE)
  X_train = X[int,]
  Y_train = Y[int]
  
  fit.cv <- cv.glmnet(X_train, Y_train, nfolds=5, alpha=0)
  
  co <- as.vector(coef(fit.cv, s = "lambda.min")) # as.vector makes an easier object to work with
  
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

ridge_bootstrap_sample(X, Y, X, Y, TRUE)$cor


ridge_boot = function(X, Y, ext.X, ext.Y, n_bootstraps=1000, method = "pearson") {
  # run many bootstraps
  # output: - vector with correlations 
  #         - matrix with features coefficients
  #         - vector with SEM
  cor_vec <- rep(NA, n_bootstraps)
  co_matrix <- matrix(NA, nrow = dim(X)[2], ncol = n_bootstraps)
  rownames(co_matrix) <- colnames(X) 
  MSE_vec <- rep(NA, n_bootstraps)
  
  for (i in c(1:n_bootstraps)) {
    out <- ridge_bootstrap_sample(X, Y, ext.X, ext.Y, TRUE, method = method)
    # browser()
    cor_vec[i] <- out$cor
    co_matrix[,i] <- out$co[-1]
    MSE_vec[i] <- out$MSE
    cat(i, "")
  } 
  return(list(cor_vec=cor_vec, co_matrix=co_matrix, MSE_vec = MSE_vec))
}

# 1000 bootstrap fits
# RUN: r_b_obj_771_p
set.seed(123)
r_b_obj_771_p <- ridge_boot(X, Y, X, Y)
head(r_b_obj_771_p$co_matrix)[,1:6]
save(r_b_obj_771_p, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_p.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_p.RData")

# RUN: r_b_obj_771_ROR_p
set.seed(123)
r_b_obj_771_ROR_p <- ridge_boot(X, Y, X, Y)
head(r_b_obj_771_ROR_p$co_matrix)[,1:6]
save(r_b_obj_771_ROR_p, file="/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_ROR_p.RData")
load("/Users/anders/Documents/MASTER/Cancer/R_codeP01/instances/r_b_obj_771_ROR_p.RData")



# rb_object
histogram(r_b_obj_771_p$cor_vec, breaks = 99,
          xlab = "Proliferation score", 
          main = "Ridge, (trail 1, arm Letro+Ribo)")





# Various objects
save(rb_object, file="rb_object_6genes.RData")
save(rb_object, file="rb_object_nodes.RData")
save(rb_object, file="rb_object_AllGenes_Prolif.RData")
save(rb_object, file="rb_object_ALLGenes_ROR.RData")

# stored object can be loaded to save time
load("rb_object_6Genes.RData")
load("rb_object_Nodes.RData")
load("rb_object_AllGenes_Prolif.RData")
load("rb_object_ALLGenes_ROR.RData")

## Summery result of lasso bootstrap object
# plot cross-validated error as a function of lambda and alpha

# Correlation
cor_vec <- rb_object$cor_vec
sum(is.na(cor_vec))/length(cor_vec)    # fraction of NA
mean(cor_vec, na.rm=TRUE)
var(cor_vec, na.rm=TRUE)
par(mfrow=c(1,1))
hist(cor_vec, breaks = 100)

# MSE
MSE_vec <- rb_object$MSE_vec
mean(MSE_vec)
var(MSE_vec)
hist(MSE_vec, breaks=100)

# Features
genes_mean <- rowMeans(rb_object$co_matrix)
apply(rb_object$co_matrix, 1, var)
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


## Analyzing selected features (genes/nodes...) in lb_object
# count the presence of the individual features for all bootstrap models
covariates_n <- 6  # If 6 genes
covariates_n <- 771  # If ALL genes
covariates_n <- 8    # If nodes
vector_1 <- c(1: covariates_n)
vector_2 <- rb_object[[2]] - 1 # shift numbers so intercept becomes 0 and first covariate is 1
covariates_count <- rowSums(outer(vector_1, vector_2, "=="))
# put gene names on the covariates_count vector
covariates_w_names <- setNames(covariates_count, colnames(X))
show(covariates_w_names)

# extract covariates selected at least once
covariates_w_names_ind <- which(covariates_w_names!=0)
covariates_w_names <- covariates_w_names[covariates_w_names_ind]
show(covariates_w_names)
max(covariates_w_names)
ind <- which(covariates_w_names == max(covariates_w_names))
covariates_w_names[ind]
covariates_w_names[order(covariates_w_names, decreasing = TRUE)]


# extract covariates selected at certain amount of times
genes_of_interest <- function(vector, times_selected, above=TRUE){
  # return: names of genes selected a distinct numbers of times
  if (max(vector) < times_selected){
    return(print("argument times_selected has to high value"))
  }
  if (above == TRUE){inds <- which(vector >= times_selected)}
  else{inds <- which(vector == times_selected)}
  
  variable_names <- names(inds)
  return(variable_names)
}



test_genes = genes_of_interest(covariates_w_names, 100, TRUE)
show(test_genes)

fm_test_genes = as.formula(paste("Proliferation.Score", "~", paste(test_genes, collapse = "+")))
lm(fm_test_genes, data = df04)

#############################################
## Bellow is code I dont use so much anymore

B <- matrix(NA,nrow = 2, ncol = 4)
colnames(B) <- colnames(X) 
# Set row names
rownames(B) <- c("Row 1", "Row 2")
rownames(B) <- paste0("Row ", 1:nrow(B)) # Equivalent


# count the number of times the selected covariates are present in the models
covariates_count_selected <- rowSums(outer(c(1:max(covariates_w_names)), covariates_w_names, "=="))
names(covariates_count_selected) <- c(1:max(covariates_w_names))
show(covariates_count_selected)



regression_bootstrap_test_genes <- function(df, formula, boot_fraction=0.632){
  N = nrow(df)
  
  int <- sample(N, size = N*boot_fraction, replace = TRUE)
  train <- df[int,]
  test <- df[-int,]
  
  lm <- lm(formula, data = train)
  pred_values <- predict(lm, newdata = test)
  out <- cor(pred_values, test$Proliferation.Score)
  return(out)
}

regression_cor_boot_test_genes = function(df, formula, boot_fraction=0.632, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] <- regression_bootstrap_test_genes(df, formula)
  }  
  return(cor_vec)
}

regression_test_genes <- regression_cor_boot_test_genes(df04, fm_test_genes, boot_fraction=0.632, 1000)

mean(regression_test_genes, na.rm=TRUE)
var(regression_test_genes, na.rm=TRUE)
par(mfrow=c(1,1))
hist(regression_test_genes, breaks = 50, xlim = c(-0.5, 0.8))



#######


