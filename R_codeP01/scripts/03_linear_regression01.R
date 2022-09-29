## Ordinary linear regression for the 6 genes as predictors
## Prediction of SUR Proliferation score is compared to measured Proliferation score at SUR

# formulas for linear model w/o and w interaction
fm01 = sur_ProliferationScore ~ scr_CCND1 + scr_CCNE1 + scr_CDKN1A + scr_ESR1 + scr_MYC + scr_RB1
fm02 = sur_ProliferationScore ~ scr_CCND1*scr_RB1 + scr_CCND1*scr_CDKN1A + scr_CCND1*scr_ESR1 + 
  scr_CCNE1*scr_CDKN1A + scr_CCNE1*scr_RB1 + scr_MYC*scr_ESR1 + scr_MYC*scr_CDKN1A 


# simple linear model
# here: correlation against a best fit of it self - overfit...:)
lm01 = lm(fm01, data = df02_LR)
summary(lm01)

pred01 = predict(lm01, df02_LR)
cor(pred01, df02_LR$sur_ProliferationScore)
plot(pred01, df02_LR$sur_ProliferationScore)
abline(lm(df02_LR$sur_ProliferationScore ~ pred01))

# linear model with interactions
lm02 = lm(fm02, data = df02_LR)
summary(lm02)

pred02 = predict(lm02, df02_LR)
cor(pred02, df02_LR$sur_ProliferationScore)
plot(pred02, df02_LR$sur_ProliferationScore)
abline(lm(df02_LR$sur_ProliferationScore ~ pred02))


# Leave one observation out for testing and going through all observations
lm_loo <- function(formula, df){
  pred_values <- NULL
  for(i in 1:nrow(df)){
    test <- df[i,]
    train <- df[-i,]
    lm01 <- lm(formula, data = train)
    pred_values[i] <- predict(lm01, newdata = test)
  }
  return(pred_values)
}

pred_values_fm01 <- lm_loo(fm01, df02_LR)
pred_values_fm02 <- lm_loo(fm02, df02_LR)

# Result of simple linear model:
par(mfrow=c(1,2))
cat("Correlation for simple linear model: ")
print(cor(pred_values_fm01, df02_LR$sur_ProliferationScore))
cat("MSE for simple linear model: ")
print(mean((pred_values_fm01 - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_fm01, df02_LR$sur_ProliferationScore, main = "Simple linear model")

# Result of linear model with interactions:
cat("Correlation for linear model with interactions: ")
print(cor(pred_values_fm02, df02_LR$sur_ProliferationScore))
cat("MSE for linear model with interactions: ")
print(mean((pred_values_fm02 - df02_LR$sur_ProliferationScore)**2))
plot(pred_values_fm02, df02_LR$sur_ProliferationScore, main = "Linear model with interaction terms")
par(mfrow=c(1,1))


# boostrap sample for train and leave out  set for testing
regression_bootstrap <- function(df, formula){
  N = nrow(df)
  
  int <- sample(N, size = N*0.632, replace = TRUE)
  train <- df[int,]
  test <- df[-int,]
  
  lm01 <- lm(formula, data = train)
  pred_values <- predict(lm01, newdata = test)
  out <- cor(pred_values, test$sur_ProliferationScore)
  return(out)
}

regression_cor_boot = function(df, formula, n_bootstraps){
  cor_vec <- rep(NA, n_bootstraps)
  for (i in c(1:n_bootstraps)) {
    cor_vec[i] <- regression_bootstrap(df, formula)
  }  
  return(cor_vec)
}

regression_cor_boot <- regression_cor_boot(df02_LR, fm01, 1000)

mean(regression_cor_boot, na.rm=TRUE)
var(regression_cor_boot, na.rm=TRUE)
hist(regression_cor_boot, breaks = 50)



#######################################################
## some not used code that might be considered later
#######################################################

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


