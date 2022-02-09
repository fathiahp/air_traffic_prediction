########################################
# Sydney Kingsford Smith Airport Traffic
########################################

#------------
# Import data
#------------

SydneyTraf <- read.csv2('SydneyAirport2020.csv')
View(SydneyTraf)
str(SydneyTraf)

traf <- ts(SydneyTraf[,3], start=c(2009,1), frequency=12)
time(traf)

#-------------------
# Stationarity issue
#-------------------

# raw data

plot.ts(traf)
acf(ts(traf, frequency=1), main="Autocorrelogram main series")
pacf(ts(traf, frequency=1), main="Partial autocorrelogram main series")

# 1st transformation: remove increasing variance effect

ltraf <- log(traf)
plot.ts(ltraf)
acf(ts(ltraf, frequency=1), main="Autocorrelogram  log(series)")
pacf(ts(ltraf, frequency=1), main="Partial autocorrelogram log(series)")

# 2nd transformation: remove trend

dltraf <- diff(ltraf,1)
plot.ts(dltraf, xlim=c(2009,2019))
acf(ts(dltraf, frequency=1), main="Autocorrelogram  1st difference of log(series)")
pacf(ts(dltraf, frequency=1), main="Partial autocorrelogram 1st difference of log(series)")

# 3rd transformation: remove seasonality

dltraf_12 <- diff(dltraf,12)
plot.ts(dltraf_12)
acf(ts(dltraf_12, frequency=1), main="Autocorrelogram  1st difference of log(series) wo seasonality")
pacf(ts(dltraf_12, frequency=1), main="Partial autocorrelogram 1st difference of log(series) wo seasonality")

#---------------------------------
# Identification of orders p and q
#---------------------------------

#------------------------------------------------------------------------
# Estimation of a multiplicative SARIMA(1,1,1)(1,1,1) with seasonality 12 
#------------------------------------------------------------------------

mod1 <- arima(ltraf, c(1,1,1), seasonal=list(order=c(1,1,1), period=12), method='ML')
mod1

#plot of the fitted value
install.packages("TSA")
library("TSA")
fit1 <- fitted(mod1)
fit1
plot.ts(cbind(ltraf,fit1), plot.type='single', col=c('black','red'))

#--------------------
# Validation of mod1
#--------------------

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Significance of coefficients: pvalue of Student test on each coeff
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Ho: coeff=0 against H1: coeff#0
# pvalue< 5%, we accept H1

mod1$coef  #value of coefficients
mod1$var.coef #variance of coefficients
tstat <- mod1$coef/sqrt(diag(mod1$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Residuals analysis: assumption of Gaussian White Noise
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''
res1 <- mod1$residuals
plot(res1)

# Autocorrelation of the residuals
#'''''''''''''''''''''''''''''''''

acf(ts(res1,frequency=1))

# Box-Pierce test ou Ljung-Box test
Box.test(res1, lag=20, type="Ljung-Box")

# Normal Distribution assumption
#'''''''''''''''''''''''''''''''

res_norm <- res1/sqrt(mod1$sigma2) # normalized residuals
summary(res_norm)


plot(res_norm)
abline(a=2,b=0,col="red")
abline(a=-2, b=0, col="red")

#QQplot 
# check whether the points are close to the line
qqnorm(res1)
qqline(res1)


#Shapiro test
shapiro.test(res1)
# pvalue < 5% so we reject the normality assumption


# Constant variance assumption
#'''''''''''''''''''''''''''''
sq.res <- (res1)^2
acf(ts(sq.res, frequency=1))

# TSA package to apply McLeod Li test
Htest <- McLeod.Li.test(mod1, plot=FALSE)
Htest



#------------
# Prediction
#-----------


# Check of the quality of the fit wrt confidence bounds
#-------------------------------------------------------

cb80 <- mod1$sigma2^.5*qnorm(0.9) # confidence bound of security 80%
plot(cbind(ltraf, fit1-cb80, fit1+cb80), plot.type='single', 
     lty=c(1,2,2), xlim=c(2009,2019))

# Proportion of points in the confidence interval
indi <- (ltraf-(fit1-cb80))>0&(fit1+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop

# if prop > 80%, then the fit is considered good


# Validation set approach
#'''''''''''''''''''''''''

#Idea: split the sample into 2 subsamples: training set and test set

data.train <- window(ltraf, start=c(2009,1), end=c(2016,12))
# 96 obs for data.train
data.test <- window(ltraf, start=c(2017,1), end=c(2019,12))
# 36 obs for data.test

mod1.train <- arima(data.train, c(1,1,1), seasonal=list(order=c(1,1,1), period=12), method='ML')
pred1.test <- predict(mod1.train, n.ahead=174)

install.packages("forecast")
library("forecast")
accuracy(pred1.test$pred,data.test)

# plot comparing observed values and prediction of the traffic
ts.plot(traf, xlim=c(2016,2020))
lines(2.718^(pred1.test$pred), col="red")
lines(2.718^(pred1.test$pred-1.96*pred1.test$se), col=4, lty=2)
lines(2.718^(pred1.test$pred+1.96*pred1.test$se), col=4, lty=2)



#--------------------------
# Estimation of a 2nd model
#--------------------------

# 1st idea: remove the non significant coeff AR1 and SAR1

# Without SAR1
mod2 <- arima(ltraf, c(1,1,1), seasonal=list(order=c(0,1,1), period=12), method='ML')
mod2
mod1

mod2$coef  #value of coefficients
mod2$var.coef #variance of coefficients
tstat <- mod2$coef/sqrt(diag(mod2$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue

# Check for residuals autocorrelation for model2
res2 <- mod2$residuals
acf(ts(res2,frequency=1)) 

# Without AR1 coefficient
mod3 <- arima(ltraf, c(0,1,1), seasonal=list(order=c(0,1,1), period=12), method='ML')
mod3

mod3$coef  #value of coefficients
mod3$var.coef #variance of coefficients
tstat <- mod3$coef/sqrt(diag(mod3$var.coef)) # test statistic of the Student test
pvalue <- 2*(1-pnorm(abs(tstat))) #pvalue of the test 
tstat
pvalue
# both MA1 and SMA1 are highly significant!

# Check for residuals autocorrelation for model3
res3 <- mod3$residuals
plot.ts(res3)
acf(ts(res3,frequency=1)) 

# Check for residuals normality

res_norm3 <- res3/sqrt(mod3$sigma2) # normalized residuals
summary(res_norm3)

# If the residuals follow a Normal distribution, the values of res_norm
# should lie in between -2 and 2, with 95% of chance

plot(res_norm3)
abline(a=2,b=0,col="red")
abline(a=-2, b=0, col="red")

# Check for residuals constant variance
sq.res3 <- (res3)^2
acf(ts(sq.res3, frequency=1))

# Check of the quality of the fit wrt confidence bounds
#-------------------------------------------------------
install.packages("TSA")
library("TSA")
fit3 <- fitted(mod3)
fit3

cb80 <- mod3$sigma2^.5*qnorm(0.9) # confidence bound of security 80%
plot(cbind(ltraf, fit3-cb80, fit3+cb80), plot.type='single', lty=c(1,2,2), xlim=c(2000,2010))

# Proportion of points in the confidence interval
indi <- (ltraf-(fit3-cb80))>0&(fit3+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop

# Prediction
#'''''''''''''''''''''''''''''''''''''''''''''''''

pred4 <- predict(mod3, n.ahead=12)
ts.plot(traf, 2.718^pred4$pred,log="y", lty=c(1,3))






