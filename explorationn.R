library(forecast)
library(tseries)
library(zoo)
library(TSA)

ajgr <- read.csv("data/co2_ajgr.csv", header=TRUE, sep=',')

d <- as.POSIXlt(ajgr$timestamp)
co2 <- zoo(ajgr$CO2, d)
plot(co2,  xlab="Date", ylab="CO2 concentration")


#-------------- Stationnarity Tests --------------#
kpss.test(co2)
adf.test(co2)

#-------------- Raw Data --------------#
par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(2,1))
plot(co2, xlab="Time", ylab="CO2 concentration"); title('a) Raw Data', line=0.3)
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
plot(diff(diff(co2, lag=24)), xlab="Time", ylab="Data diff"); title('a) Diff Data', line=0.3)
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
best_model_parameters <- NULL
best_model_fit <- NULL
best_model_predictions <- NULL

for (s in sarimas){
  model_fit <- Arima(co2[1:576], order=s$order, seasonal=s$seasonal)
  predictions <- forecast(model_fit, h=168, level=95)[["mean"]]

  avg_error <- mean((true_values - predictions)^2)
  cat("(", s$order, ") x (", s$seasonal$order, ")_", s$seasonal$period, "\n")
  cat("loss =", avg_error, "\n")
  cat("aic =", model_fit$aic, "\n")
  cat("bic =", model_fit$bic, "\n")
  cat("loglik =", model_fit$loglik, "\n")
  
  if(is.null(best_error) || avg_error < best_error){
    best_model_parameters = s
    best_error <- avg_error
    best_model_fit = model_fit
    best_model_predictions <- predictions
  }
}

# let's show how good we predicted the last week compared to true values

plot(true_values, type='l', col='blue', xlab="Date", ylab="CO2 concentration")
title('true values vs predictions of last week')
lines(best_model_predictions, col='red')
legend(d[577:744][135], 354, legend=c("True", "Predicted"), col=c("blue", "red"), lty=1:1, cex=0.7)


# let's analyze the best model

best_full_model_fit <- Arima(as.ts(co2), order=best_model_parameters$order, seasonal=best_model_parameters$seasonal)
e <- best_full_model_fit$residuals
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
future_values <- forecast(best_full_model_fit, h=120, level=95)

from <- as.double(d[1])
from_pred <- as.double(tail(d, n=1))
to <- from_pred + 60 * 60 * 119

time <- seq(from_pred, to, by=60*60)

plot(
  zoo(future_values$fitted, d), 
  xlim=c(from, to), 
  ylim=c(min(min(co2), min(future_values$lower)), max(max(co2), max(future_values$upper))),
  xlab="Date", 
  ylab="CO2 concentration"
)
polygon(c(time, rev(time)), c(future_values$upper, rev(future_values$lower)), col = "grey", border = NA)
lines(future_values$mean, col="blue")
  
