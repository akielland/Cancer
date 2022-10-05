
# first: correlation against a best fit of it self - overfit...:)

# linear base learner
boost_m01 = glmboost(fm01, data = df02_LR)
coef(boost_m01, which = "")
par(mfrow=c(1,2))
plot(boost_m01, off2int=TRUE)
plot(boost_m01)

pred_boost01 = predict(boost_m01, df02_LR)

cor(pred_boost01, df02_LR$sur_ProliferationScore)
plot(pred_boost01, df02_LR$sur_ProliferationScore)
abline(lm(df02_LR$sur_ProliferationScore ~ pred_boost01))


# smooth - P-spline as base learner
spline01 <- gamboost(fm01, data = df02_LR,
                     baselearner = "bbs", # dfbase=4
                     control = boost_control(mstop = 50))

coef(spline01, which = "")
par(mfrow=c(2,3))
plot(spline01, off2int=TRUE)

pred_spline01 = predict(spline01, df02_LR)

cor(pred_spline01, df02_LR$sur_ProliferationScore)
par(mfrow=c(1,1))
plot(pred_spline01, df02_LR$sur_ProliferationScore)
abline(lm(df02_LR$sur_ProliferationScore ~ pred_spline01))

# boosting with stumps
stump01 = gamboost(fm01, data = df02_LR, 
                   family=Gaussian(), 
                   baselearner='btree',
                   boost_control(mstop = 100))

cols_chosen = unique(stump01$xselect())
colnames(df02)[cols_chosen+2]

par(mfrow=c(2,3))
plot(stump01)


pred_boost01 = predict(boost_m01, df02_LR)

cor(pred_boost01, df02_LR$sur_ProliferationScore)
plot(pred_boost01, df02_LR$sur_ProliferationScore)
abline(lm(df02_LR$sur_ProliferationScore ~ pred_boost01))
