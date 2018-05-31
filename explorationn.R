library(forecast)
library(tseries)
library(zoo)
library(TSA)

ajgr <- read.csv("data/co2_ajgr.csv", header=TRUE, sep=',')

d <- as.POSIXlt(ajgr$timestamp)
z <- zoo(ajgr$CO2, d)

plot(z,  xlab="Date", ylab="CO2 concentration")

#-------------- Raw Data --------------#
co2 <- as.ts(z)
kpss.test(co2)

par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(3,1))
plot(z, xlab="Time", ylab="CO2 concentration"); title('a) Raw Data', line=0.3)
Acf(co2, lag.max = 100); title('b) ACF', line=0.3)
spectrum(ajgr$CO2, ylab="Periodogram", xlab="Frequency", main='')


#-------------- Log data --------------#
co2_log <- log(co2-300)
plot(co2_log)
kpss.test(co2_log)

#-------------- Diff Data --------------#
co2_d <- diff(co2, lag=24)
kpss.test(co2_d)

#-------------- Diff^2 data --------------#
co2_dd <- diff(diff(co2, lag=24))
kpss.test(co2_dd)

par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(3,1))
plot(diff(diff(z, lag=24)), xlab="Time", ylab="Data diff"); title('a) Diff Data', line=0.3)
Acf(co2_dd, lag.max = 100); title('b) ACF', line=0.3)
Pacf(co2_dd, lag.max = 100); title('c) PACF', line=0.3)


# OUR MODEL
co2_choosen <- co2_dd

Acf(co2_choosen, lag.max = 100)
Pacf(co2_choosen, lag.max = 100)

# => SAR
fit0 <- Arima(co2, order=c(2, 1, 2), seasonal=list(order=c(0, 1, 1), period=24))
e <- fit0$residuals
plot(e)
Acf(e)
Pacf(e)
cpgram(e)
qqnorm(e) #not good because we have extrem values

# Ljung-Box
lbt <- c(); for (h in 6:25) lbt[h] <- Box.test(e,lag=h,type='Ljung-Box',fitdf=5)$p.value
plot(lbt, ylim=c(0,1)); abline(h=0.05,col='blue',lty='dotted')

# forecasting
fit0 <- Arima(co2[0:620], order=c(2, 1, 2), seasonal=list(order=c(0, 1, 1), period=24))
plot(forecast(fit0, h=500, level=95))

fit1 <- Arima(co2[0:576], order=c(2, 1, 2), seasonal=list(order=c(0, 1, 1), period=24))
f <- forecast(fit1, h=168, level=95)

valuesForcasted <- f[["mean"]]
trueValues <- co2[577:744]

x <-  seq(1, 168, 1)

plot(x, trueValues, type='l', col='blue', ylab="CO2 concentration")
title('Predicted values vs true values')
legend(135, 354, legend=c("True", "Predicted"), col=c("blue", "red"),
       lty=1:1, cex=0.7)
#axis(1, at=(1:length(x):10), labels = d[577:744])
lines(x, valuesForcasted, col='red')