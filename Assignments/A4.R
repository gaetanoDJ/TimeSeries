library(FinTS)
library(tseries)
library(aTSA)

attach(irates)
ts(r1)
ts(r3)
tseries::adf.test(r1)
tseries::adf.test(r3)
### r1 and r3 contain a unit root

Reg1 <-lm(r1 ~ r3)
Reg2 <-lm(r3 ~ r1)

summary(Reg1)
summary(Reg2)
coint.test(r1, r3)
coint.test(r3,r1)

# pvalue greater than 0.05 so we reject the null hypothesis that the two time series are not cointegrated


bfx <- as.matrix(cbind(r1,r3), demean=FALSE)
po.test(bfx)

bfx1 <- as.matrix(cbind(r3,r1), demean=FALSE)
po.test(bfx1)


# po.test also states that we fail to reject the null hypothesis that there is no cointegration.
