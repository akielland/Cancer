library(MASS)
library(tidyverse)
# install.packages("tidyverse")

# Importing data
df01 <- readxl::read_xlsx('/Users/anders/Documents/Master/data/SCR_SUR_6genes_noNAN.xlsx')

df01 <- read.delim2('/Users/anders/Documents/MASTER/Cancer/R_codeP01/data/pro')

# Make dataframe with only timepoint SCR (both interventions)
df01.SCR <- df01[df01$timepoint %in% "SCR", ]

# Split in training and test set
n.obs = nrow(df01.SCR)
set.seed(1)
sample = sample(n.obs, 60)
df01.SCR_train = df01.SCR[sample,]
df01.SCR_test = df01.SCR[-sample,]

formula01 = ProliferationScore ~ CCND1 + CCNE1 + CDKN1A + ESR1 + MYC + RB1

lm01 = lm(formula01, data = df01.SCR)
summary(lm01)

lm02 = lm(formula01, data = df01.SCR_train)
summary(lm02)

pred01 = predict(lm02, df01.SCR_test)
cor(pred01, df01.SCR_test$ProliferationScore)
plot(pred01, df01.SCR_test$ProliferationScore)
abline(lm(df01.SCR_test$ProliferationScore ~ pred01))
cor(pred01[-17], df01.SCR_test$ProliferationScore[-17])
plot(pred01, df01.SCR_test$ProliferationScore, xlim = c(0, 0.6))
abline(lm(df01.SCR_test$ProliferationScore[-17] ~ pred01[-17]))


formula02 = ProliferationScore ~ CCND1*RB1 + CCND1*CDKN1A + CCND1*ESR1 + CCNE1*CDKN1A +  CCNE1*RB1 + MYC*ESR1 + MYC*CDKN1A  
lm03 = lm(formula02, data = df01.SCR)
summary(lm03)
lm04 = lm(formula02, data = df01.SCR_train)
summary(lm04)

stepAIC(lm01, direction="both")
stepAIC(lm02, direction="both")
stepAIC(lm03, direction="both")
stepAIC(lm04, direction="both", trace=0, k = log(95))


df01.SCR.X <- as.matrix(subset(df01.SCR, select=-c(UniqueID, timepoint, TrialArmNeo, ProliferationScore)))
df01.SCR.Y <- as.matrix(df01.SCR["ProliferationScore"])

Lasso01.cv <- cv.glmnet(df01.SCR.X, df01.SCR.Y)
plot(Lasso01.cv)
coef(Lasso01.cv, s="lambda.min")
