library(forecast)
library(tseries)
library(zoo)
library(TSA)

ajgr <- read.csv("data/co2_ajgr.csv", header=TRUE, sep=',')

d <- as.POSIXlt(ajgr$timestamp)
z <- zoo(ajgr$CO2, d)

plot(z,  xlab="Date", ylab="CO2 concentration")

#-------------- Raw Data --------------#
co2 <- ajgr$CO2
kpss.test(co2)
adf.test(co2)

par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(2,1))
plot(z, xlab="Time", ylab="CO2 concentration"); title('a) Raw Data', line=0.3)
Acf(co2, lag.max = 100); title('b) ACF', line=0.3)

spectrum(ajgr$CO2, ylab="Periodogram", xlab="Frequency", main='')
title('c) Spectrum', line=0.3)


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
adf.test(co2_dd)

par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(3,1))
plot(diff(diff(z, lag=24)), xlab="Time", ylab="Data diff"); title('a) Diff Data', line=0.3)
Acf(co2_dd, lag.max = 100); title('b) ACF', line=0.3)
Pacf(co2_dd, lag.max = 100); title('c) PACF', line=0.3)


#-------------- OUR MODEL --------------#
co2_choosen <- co2_dd

Acf(co2_choosen, lag.max = 100)
Pacf(co2_choosen, lag.max = 100)

# many sarimas models

sarimas <- list(
  list(order=c(2, 1, 2), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(1, 1, 2), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(2, 1, 1), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(1, 1, 1), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(1, 1, 0), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(0, 1, 1), seasonal=list(order=c(0, 1, 1), period=24))
)


#-------------- Comparison of models --------------#
# let's find the best model by predicting the last week
# of data we have using the weeks before

true_values <- co2[577:744]

best_error <- NULL
best_model_fit <- NULL
best_model_predictions <- NULL

for (s in sarimas){
  model_fit <- Arima(co2[1:576], order=s$order, seasonal=s$seasonal)
  predictions <- forecast(model_fit, h=168, level=95)[["mean"]]

  avg_error <- mean((true_values - predictions)^2)
  print(avg_error)
  
  if(is.null(best_error) || avg_error < best_error){
    print("new best model changed")
    print(paste("it has an parameters =", s))
    print(paste("it has an average error = ", avg_error))
    best_error <- avg_error
    best_model_fit = model_fit
    best_model_predictions <- predictions
  }
}


# let's show how good we predicted the last week compared to true values

x <-  seq(1, 168, 1)

plot(x, true_values, type='l', col='blue', ylab="CO2 concentration")
title('Predicted values vs true values')
lines(x, best_model_predictions, col='red')
legend(135, 354, legend=c("True", "Predicted"), col=c("blue", "red"), lty=1:1, cex=0.7)



# let's analyze the best model

e <- best_model_fit$residuals
plot(e)
Acf(e)
Pacf(e)
cpgram(e)
qqnorm(e)
qqline(e, col='red') #not good because we have extrem values


#-------------- Ljung-Box --------------#
lbt <- c(); for (h in 6:25) lbt[h] <- Box.test(e,lag=h,type='Ljung-Box',fitdf=5)$p.value
plot(lbt, ylim=c(0,1)); abline(h=0.05,col='blue',lty='dotted')



#-------------- Forecasting future value --------------#
plot(forecast(best_model_fit, h=120, level=95))


