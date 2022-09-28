## Ordinary linear regression for the 6 genes as predictors
## Prediction of SUR Proliferation score is compared to measured Proliferation score at SUR

# formulas for linear model w/o and w interaction
fm01 = sur_ProliferationScore ~ scr_CCND1 + scr_CCNE1 + scr_CDKN1A + scr_ESR1 + scr_MYC + scr_RB1
fm02 = sur_ProliferationScore ~ scr_CCND1*scr_RB1 + scr_CCND1*scr_CDKN1A + scr_CCND1*scr_ESR1 + 
  scr_CCNE1*scr_CDKN1A + scr_CCNE1*scr_RB1 + scr_MYC*scr_ESR1 + scr_MYC*scr_CDKN1A 

# creating train and test set
n_holdout_test = 0
int <- sample.int(nrow(df02_LR), size = nrow(df02_LR) - n_holdout_test, replace = FALSE)
train = df02_LR[int,]
test = df02_LR[-int,]

# simple linear model
lm01 = lm(fm01, data = train)
summary(lm01)

pred01 = predict(lm01, test)
cor(pred01, test$sur_ProliferationScore)
plot(pred01, test$sur_ProliferationScore)
abline(lm(test$sur_ProliferationScore ~ pred01))

# linear model with interactions
lm02 = lm(fm02, data = train)
summary(lm02)

pred02 = predict(lm02, test)
cor(pred02, test$sur_ProliferationScore)
plot(pred02, test$sur_ProliferationScore)
abline(lm(test$sur_ProliferationScore ~ pred02))


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


