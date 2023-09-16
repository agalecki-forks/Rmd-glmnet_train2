names(T.1)

cvfitT1 <- T.1$cvfit
myglance(cvfitT1)
mytidy(cvfitT1)
xnewT1 <- T.1$xnew
dim(xnewT1)

ySurvT1 <- T.1$ySurv
dim(ySurvT1)

fitT3 <- T.3$fit

predT3 <- predict(fitT3, newx = xnewT1)


