library(forecast)
library(tseries)
library(zoo)
library(TSA)


#-------------- Raw Data Loading --------------#
ajgr <- read.csv("data/co2_ajgr.csv", header=TRUE, sep=',')

d <- as.POSIXlt(ajgr$timestamp)
co2 <- zoo(ajgr$CO2, d)
plot(co2,  xlab="Time", ylab="CO2 concentration")


#-------------- Stationnarity Tests --------------#
kpss.test(co2)
adf.test(co2)

#-------------- Raw Data --------------#
par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(2,1))
plot(co2, xlab="Time", ylab="CO2 ppm"); title('Raw Data', line=0.3)
Acf(co2, lag.max = 100); title('ACF of raw data', line=0.3)

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
plot(diff(diff(co2, lag=24)), xlab="Time", ylab="Data diff"); title('Data differentiated', line=0.3)
Acf(co2_dd, lag.max = 100); title('ACF of data differentiated', line=0.3)
Pacf(co2_dd, lag.max = 100); title('PACF of data differentiated', line=0.3)


#-------------- Stationnary model we will fit --------------#
co2_choosen <- co2_dd

Acf(co2_choosen, lag.max = 100)
Pacf(co2_choosen, lag.max = 100)


#-------------- Comparison of models --------------#
# let's find the best model by predicting the last week
# of data we have using the weeks before

# many possible sarimas models
sarimas <- list(
  list(order=c(2, 1, 2), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(1, 1, 2), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(2, 1, 1), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(1, 1, 1), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(1, 1, 0), seasonal=list(order=c(0, 1, 1), period=24)),
  list(order=c(0, 1, 1), seasonal=list(order=c(0, 1, 1), period=24))
)

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

par(mfrow=c(1,1))
plot(true_values, type='l', col='blue', xlab="Time", ylab="CO2 ppm")
title('true values vs predictions of last week')
lines(best_model_predictions, col='red')
legend(d[577:744][145], 354, legend=c("True", "Predictions"), col=c("blue", "red"), lty=1:1, cex=0.7)


#-------------- Final Model Analysis --------------#



par(mar=c(3,3,1.2,0.1), mgp=c(1.5,0.4,0))
par(mfrow=c(3,1))

best_full_model_fit <- Arima(as.ts(co2), order=best_model_parameters$order, seasonal=best_model_parameters$seasonal)
e <- zoo(best_full_model_fit$residuals, d)
plot(e, xlab="Time", ylab="Residuals"); title('Residuals after model fitting')
Acf(e, lag.max = 100); title("ACF of residuals")
Pacf(e, lag.max = 100); title("PACF of residuals")

#-------------- Ljung-Box --------------#
par(mfrow=c(1,1))

lbt <- c(); for (h in 6:25) lbt[h] <- Box.test(e,lag=h,type='Ljung-Box',fitdf=5)$p.value
plot(lbt, ylim=c(0,1)); abline(h=0.05,col='blue',lty='dotted'); title("Ljung-Box Test of residuals")


par(mfrow=c(1,3))
cpgram(e, main="Cum. Periodogram of residuals")

qq <- qqnorm(e); qqline(e, col='red') #not good because we have extrem values
all <- qq$y
qq$y[[which.min(qq$y)]]
max(all)
range <- c(qq$x[[which.min(qq$x)]], qq$x[[which.max(qq$x)]])
plot(qq$x, qq$y, xlim=range, ylim=range, xlab="Theoretical Quantiles", ylab="Sample Quantiles"); title("Normal Q-Q Plot with th. qu. ranges"); qqline(e, col='red')


#-------------- Forecasting future value --------------#

par(mfrow=c(1,1))
pred_length = 24 * 7
future_values <- forecast(best_full_model_fit, h=pred_length, level=95)

from <- as.double(d[1])
from_pred <- as.double(tail(d, n=1))
to <- from_pred + 60 * 60 * (pred_length - 1)

time <- seq(from_pred, to, by=60*60)

plot(
  zoo(future_values$fitted, d), 
  xlim=c(from, to), 
  ylim=c(min(min(co2), min(future_values$lower)), max(max(co2), max(future_values$upper))),
  xlab="Time", 
  ylab="CO2 concentration"
)
title("CO2 concentration predictions for one week")
polygon(c(time, rev(time)), c(future_values$upper, rev(future_values$lower)), col = "grey", border = NA)
lines(future_values$mean, col="blue")
  
