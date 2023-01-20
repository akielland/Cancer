################
## Post Lasso ##
################

## Lasso is used for feature selection
## Test the features with:
## - Linear regression
## - Ridge regression

test_genes01 <- features_ordered.topp_list[1]$rowname[1:30]
length(test_genes01)
test_genes01

fm_test01 = as.formula(paste("Y", "~", paste(test_genes01, collapse = "+")))
fm_test01

fit.lm <- lm(fm_test01, data = Proliferation_ALLgenes)
fit.lm
ext.pred_ROR.lm <- predict(fit.lm, newdata = df04)

summary(ext.pred_ROR.lm)
hist(ext.pred_ROR.lm)

## Correlations
plot(ext.pred_ROR.lm, unlist(ext.PFS_months))
# Pearson
cor(ext.pred_ROR.lm, ext.PFS_months, method = "pearson")
summary(ext.PFS_months)
cor(ext.pred_ROR, ext.PFS_status)
summary(ext.PFS_status)

# Spearman
cor(ext.pred_ROR.lm, ext.PFS_months, method = "spearman")
# p-values
cor.test(ext.pred_ROR.lm, unlist(ext.PFS_months), method = "spearman")



fit.ridge <- cv.glmnet(X_train, Y_train, nfolds=5, alpha=0)






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

test_genes = genes_of_interest(covariates_w_names, 200, TRUE)
show(test_genes)

test_Y <- "Proliferation.Score"
test_Y <- Y # "ROR_P_Subtype_Proliferation"
fm_test = as.formula(paste("test_Y", "~", paste(test_genes, collapse = "+")))
fit.lm <- lm(fm_test, data = Proliferation_ALLgenes)
predict(fit.lm, data = df04)


features_ordered.topp_list



        