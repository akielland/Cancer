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

n_holdout = 0
int <- sample.int(nrow(df02_LR), size = nrow(df02_LR) - n_holdout, replace = FALSE)
train = df02_LR[int,]
test = df02_LR[-int,]

# formulas for linear model wo and w interaction
fm01 = sur_ProliferationScore ~ scr_CCND1 + scr_CCNE1 + scr_CDKN1A + scr_ESR1 + scr_MYC + scr_RB1
fm02 = sur_ProliferationScore ~ scr_CCND1*scr_RB1 + scr_CCND1*scr_CDKN1A + scr_CCND1*scr_ESR1 + 
  scr_CCNE1*scr_CDKN1A + scr_CCNE1*scr_RB1 + scr_MYC*scr_ESR1 + scr_MYC*scr_CDKN1A 

# simple linear model
lm01 = lm(fm01, data = train)
summary(lm01)

pred01 = predict(lm01, test)
cor(pred01, test$sur_ProliferationScore)
plot(pred01, test$sur_ProliferationScore)
abline(lm(test$sur_ProliferationScore ~ pred01))

lm02 = lm(fm02, data = train)
summary(lm02)

pred02 = predict(lm02, test)
cor(pred02, test$sur_ProliferationScore)
plot(pred02, test$sur_ProliferationScore)
abline(lm(test$sur_ProliferationScore ~ pred02))

#

lm_loo <- function(formula){
  pred_values <- NULL
  for(i in 1:nrow(df02_LR)){
    test <- df02_LR[i,]
    train <- df02_LR[-i,]
    lm01 <- lm(formula, data = train)
    pred_values[i] <- predict(lm01, newdata = test)
  }
  cat("Correlation: ")
  print(cor(pred_values, df02_LR$sur_ProliferationScore))
  cat("MSE: ")
  print(mean((pred_values - df02_LR$sur_ProliferationScore)**2))
  return(pred_values)
}

pred_values_fm01 <- lm_loo(fm01)
pred_values_fm02 <- lm_loo(fm02)

plot(pred_values_fm01, df02_LR$sur_ProliferationScore)
plot(pred_values_fm02, df02_LR$sur_ProliferationScore)


#######################################################
stepAIC(lm01, direction="both")
stepAIC(lm02, direction="both")
stepAIC(lm01, direction="both", trace=0, k = log(95))
stepAIC(lm02, direction="both", trace=0, k = log(95))


