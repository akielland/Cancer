# scatter plot of 6 genes at SCR vs Proliferation score at SUR

X <- dplyr::select(df02_LR,  scr_CCND1, scr_CCNE1, scr_CDKN1A, scr_ESR1, scr_MYC, scr_RB1)
Y <- select(df02_LR, sur_ProliferationScore)

par(mfrow=c(2,3))

for(i in 1:(ncol(X))) {
  x = as.matrix(X[,i])
  y = as.matrix(Y)
  plot(x, y, 
       xlab=colnames(X)[i],
       ylab=colnames(Y)[1],
       pch=20, col='blue',
       main=paste('Plot number',i))
  abline(lm(y ~ x), col = "red", lwd=3)
  mtext(paste("Correlation:", round(cor(x, y), 3)), side=1, line=-1.5, cex=1.1)
  # cat(names[i])
  # print(cor(x,y), method = c("pearson"))
}
par(mfrow=c(1,1))
