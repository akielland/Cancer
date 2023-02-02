###################################
## Ordinary linear model
## AIC
## Bootstrap
###################################
##
## Predictors: genes (various structures), node values
## Output: correlations

## FORMULAS ##
##############
# 6 genes:
fm00 <- CCND1 ~ CCND1 + CCNE1
fm01 <- Y ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1
fm02 <- Y ~ CCND1*RB1 + CCND1*CDKN1A + CCND1*ESR1 + CCNE1*CDKN1A + CCNE1*RB1 + MYC*ESR1 + MYC*CDKN1A
fm03 <- Y ~ CDKN1A + ESR1 + MYC + ESR1:MYC + CDKN1A:MYC
# nodes:
fm04 <- Y ~ cyclinD1 + cyclinD1Palbo + p21 + cyclinD1p21 + cMyc + cyclinEp21 + Rb1 + ppRb1
# 771 genes:
fm05 <- Y ~ page(genes, collapse("+"))
fm05 <- as.formula(paste(Y, "~", genes))
fm05 <- formula(paste(x, collapse = " ")) 
fm05 <- as.formula(paste("Y", paste(genes, collapse="+"), sep=" ~ "))


##############################################
## Linear models evaluated at training data ##
##############################################

# here: correlation against a best fit of it self - over-fit...:)
par(mfrow=c(1,3))

# linear models with 6 gene dataframe
lm01 = lm(fm01, data = Proliferation_6genes)
summary(lm01)
pred01 = predict(lm01, Proliferation_6genes)
cor(pred01, Proliferation_6genes$Y)
plot(pred01, Proliferation_6genes$Y)
abline(lm(Proliferation_6genes$Y ~ pred01))

# linear model with Node dataframe
lm02 = lm(fm04, data = Nodes_Proliferation)
summary(lm02)
pred02 = predict(lm02, Nodes_Proliferation)
cor(pred02, Nodes_Proliferation$Y)
plot(pred02, Nodes_Proliferation$Y)
abline(lm(Nodes_Proliferation$Y ~ pred02))


# stepwise to browse important of genes and 
# to create a minimal interaction model
stepAIC(lm01, direction="both", trace=0)
stepAIC(lm02, direction="both", trace=0)
stepAIC(lm01, direction="both", trace=0, k = log(95))
stepAIC(lm02, direction="both", trace=0, k = log(95))


# linear model with AIC reduced predictors ( with interactions)
lm03 = lm(fm03, data = df02_LR)
summary(lm03)
pred03 = predict(lm03, df02_LR)
cor(pred03, df02_LR$sur_ProliferationScore)
plot(pred03, df02_LR$sur_ProliferationScore)
abline(lm(df02_LR$sur_ProliferationScore ~ pred03))


###############################################################
# bootstrap sample for train and original sample for testing ##
###############################################################

regression_bootstrap <- function(df, formula){
  N = nrow(df)
  
  int <- sample(1:N, size = N, replace = TRUE)
  train <- df[int,]
  test <- df
  
  fit <- lm(formula, data = train)
  pred_values <- predict(fit, newdata = test)
  out <- cor(pred_values, test$Y)
  return(out)
}

regression_cor_boot = function(df, formula, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] <- regression_bootstrap(df, formula)
  }  
  return(cor_vec)
}

par(mfrow=c(3,1))

# 6 genes
regression_cor_boot_fm01 <- regression_cor_boot(Proliferation_6genes, fm01, 1000)
mean(regression_cor_boot_fm01, na.rm=TRUE)
var(regression_cor_boot_fm01, na.rm=TRUE)
hist(regression_cor_boot_fm01, breaks = 50, xlim = c(-0.5, 0.8))

regression_cor_boot_fm02 <- regression_cor_boot(Proliferation_6genes, fm02, 1000)
mean(regression_cor_boot_fm02, na.rm=TRUE)
var(regression_cor_boot_fm02, na.rm=TRUE)
hist(regression_cor_boot_fm02, breaks = 50, xlim = c(-0.5, 0.8))

regression_cor_boot_fm03 <- regression_cor_boot(Proliferation_6genes, fm03, 1000)
mean(regression_cor_boot_fm03, na.rm=TRUE)
var(regression_cor_boot_fm03, na.rm=TRUE)
hist(regression_cor_boot_fm03, breaks = 50, xlim = c(-0.5, 0.8))

# Nodes
regression_cor_boot_fm04 <- regression_cor_boot(Nodes_Proliferation, fm04, 1000)
mean(regression_cor_boot_fm04, na.rm=TRUE)
var(regression_cor_boot_fm04, na.rm=TRUE)
hist(regression_cor_boot_fm04, breaks = 50, xlim = c(-0.5, 0.8))



#####################################################################
## Leave one observation out for testing and use rest for modeling ##
#####################################################################
# (going through all observations
lm_loo <- function(formula, df){
  pred_values <- NULL
  for(i in 1:nrow(df)){
    test <- df[i,]
    train <- df[-i,]
    lm_ <- lm(formula, data = train)
    pred_values[i] <- predict(lm_, newdata = test)
  }
  return(pred_values)
}

pred_values_fm01 <- lm_loo(fm01, df02_LR)
pred_values_fm02 <- lm_loo(fm02, df02_LR)
pred_values_fm03 <- lm_loo(fm03, df02_LR)

# Result of simple linear model:
par(mfrow=c(1,2))
cat("Correlation for simple linear model (fm01): ")
print(cor(pred_values_fm01, df02_LR$sur_ProliferationScore))
cat("MSE for simple linear model: ")
print(mean((pred_values_fm01 - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_fm01, df02_LR$sur_ProliferationScore, main = "Simple linear model")

# Result of linear model with many interactions (fm02):
cat("Correlation for linear model with interactions: ")
print(cor(pred_values_fm02, df02_LR$sur_ProliferationScore))
cat("MSE for linear model with interactions: ")
print(mean((pred_values_fm02 - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_fm02, df02_LR$sur_ProliferationScore, main = "Linear model with interaction terms")
par(mfrow=c(1,1))

cat("Correlation for linear model with selected interactions (fm03): ")
print(cor(pred_values_fm03, df02_LR$sur_ProliferationScore))
cat("MSE for linear model with few interactions: ")
print(mean((pred_values_fm03 - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_fm03, df02_LR$sur_ProliferationScore, main = "Linear model with few interaction terms")




#############################################################
## boostrap sample for train and leave out set for testing ##
#############################################################

regression_bootstrap <- function(df, formula, boot_fraction=0.632){
  N = nrow(df)
  
  int <- sample(N, size = N*boot_fraction, replace = TRUE)
  train <- df[int,]
  test <- df[-int,]
  
  lm01 <- lm(formula, data = train)
  pred_values <- predict(lm01, newdata = test)
  out <- cor(pred_values, test$sur_ProliferationScore)
  return(out)
}

regression_cor_boot = function(df, formula, boot_fraction=0.632, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] <- regression_bootstrap(df, formula)
  }  
  return(cor_vec)
}

regression_cor_boot_fm01 <- regression_cor_boot(df02_LR, fm01, boot_fraction=0.632, 1000)

par(mfrow=c(3,1))
mean(regression_cor_boot_fm01, na.rm=TRUE)
var(regression_cor_boot_fm01, na.rm=TRUE)
hist(regression_cor_boot_fm01, breaks = 50, xlim = c(-0.5, 0.8))

regression_cor_boot_fm02 <- regression_cor_boot(df02_LR, fm02, boot_fraction=0.632, 1000)

mean(regression_cor_boot_fm02, na.rm=TRUE)
var(regression_cor_boot_fm02, na.rm=TRUE)
hist(regression_cor_boot_fm02, breaks = 50, xlim = c(-0.5, 0.8))

regression_cor_boot_fm03 <- regression_cor_boot(df02_LR, fm03, boot_fraction=0.632, 1000)

mean(regression_cor_boot_fm03, na.rm=TRUE)
var(regression_cor_boot_fm03, na.rm=TRUE)
hist(regression_cor_boot_fm03, breaks = 50, xlim = c(-0.5, 0.8))

# df.hist <- as.data.frame(rbind(regression_cor_boot_fm03,regression_cor_boot_fm02, regression_cor_boot_fm01))
# 
# ggplot(df.hist, aes(x=df.hist, fill=)) + 
#   geom_histogram(alpha=0.2, position="identity")
# 
# ggplot(df.hist, aes(x=df.hist)) + 
#   geom_histogram(data = df.hist[1,], fill = "red", alpha = 0.2) + 
#   geom_histogram(data = df.hist[1,], fill = "blue", alpha = 0.2) +
#   geom_histogram(data = df.hist[1,], fill = "green", alpha = 0.2) 
# 
# 
# 
# plot_multi_histogram <- function(df, feature, label_column) {
#   plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
#     geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
#     geom_density(alpha=0.7) +
# 
#   plt + guides(fill=guide_legend(title=label_column))
# }
# 
# plot_multi_histogram(df.hist, 'Sepal.Width')

#######################################################
## some old code that might be needed
#######################################################

## Ordinary linear regression for the 6 genes as predictors
## Prediction of SUR Proliferation score is compared to measured Proliferation score at SUR

# formulas for linear model w/o and w interaction
# fm01 <- sur_ProliferationScore ~ scr_CCND1 + scr_CCNE1 + scr_CDKN1A + scr_ESR1 + scr_MYC + scr_RB1
# 
# fm02 <- sur_ProliferationScore ~ scr_CCND1*scr_RB1 + scr_CCND1*scr_CDKN1A + scr_CCND1*scr_ESR1 + 
#   scr_CCNE1*scr_CDKN1A + scr_CCNE1*scr_RB1 + scr_MYC*scr_ESR1 + scr_MYC*scr_CDKN1A 
# 
# fm03 <- sur_ProliferationScore ~ scr_CDKN1A + scr_ESR1 + scr_MYC + 
#   scr_ESR1:scr_MYC + scr_CDKN1A:scr_MYC 


# stepAIC(lm01, direction="both")
# stepAIC(lm02, direction="both")
# stepAIC(lm01, direction="both", trace=0, k = log(95))
# stepAIC(lm02, direction="both", trace=0, k = log(95))

# Split in training and test set
# set.seed(001)
# iris_split <- initial_split(, prop = 0.7)
# train_data <- training(iris_split )
# test_data <- testing(iris_split )
# 
# # create an ID
# df <- df %>% mutate(id = row_number())
# train <- df %>% sample_frac(.70)
# test  <- anti_join(df, train, by = 'id')  


