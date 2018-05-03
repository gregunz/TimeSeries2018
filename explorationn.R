library(forecast)
library(tseries)

ajgr <- read.csv("data/co2_ajgr.csv", header=TRUE, sep=',')
co2 <- ts(ajgr['CO2'], frequency = 24, start = as.Date('2017-10-01 00:00:00'), end = as.Date('2017-10-31 23:00:00'))

plot(co2)

co2_d <- diff(co2)
co2_2d <- diff(diff(co2, lag=24))
co2_log <- log(co2 - 300)

plot(co2_d)
plot(co2_log)

# OUR MODEL
co2_choosen <- co2_2d

Acf(co2_choosen, lag.max = 100)
Pacf(co2_choosen, lag.max = 100)

# => SAR
fit0 <- Arima(co2_choosen, order=c(0, 1, 1), seasonal=list(order=c(0, 1, 1), period=24))
e <- fit0$residuals
plot(e)
Acf(e)
Pacf(e)
cpgram(e)

# Ljung-Box
lbt <- c(); for (h in 3:25) lbt[h] <- Box.test(e,lag=h,type='Ljung-Box',fitdf=2)$p.value
plot(lbt, ylim=c(0,1)); abline(h=0.05,col='blue',lty='dotted')

#par(mar=c(3,3,1.2,0.1),mgp=c(1.5,0.4,0))
#par(mfrow=c(3,1))
#plot(e,main='',ylab='',xlab='Time'); title('Residuals', line=0.3)
#Acf(e,main=''); title('ACF of residuals', line=0.3)
#plot(3:25, lbt[3:25], ylim=c(0,1),xlab='DF', ylab='p value',main=''); abline(h=0.05,col='blue',lty='dotted'); title('p values for Ljung-Box statistic (adjusted DF)', line=0.3)
