##################################
#                                #
## Advance Time Series Analysis ##
#                                #
# Assignment                     #
# Name: Jerko Radan Cruz         #
# r00773435                      #
#                                #
##################################

rm(list=ls())
# Import data

#Data of Microsoft stock price
mydata1<-read.table(file="MSFT.csv", sep = ",", header=TRUE,dec=".")
names(mydata1)
attach(mydata1)
head(mydata1)

full_stock_ts <- ts(Close, frequency = 12, start = c(1986, 3))

#Data of 10-year breakeven inflation rate
mydata2<-read.table(file="T10YIEM.csv", sep = ",", header=TRUE,dec=".")
names(mydata2)
attach(mydata2)
head(mydata2)

infl_expect_ts <- ts(T10YIEM, frequency = 12, start = c(2003,1))

# Accotate the time series for MSFT stocks from December 2002 to October 2019
# Reasons: match the time period from the 10-year breakeven inflation rate (+1 since it will go in differences)
stock_ts <- window(full_stock_ts, start = c(2002,12), end = c(2019,10))


###################################################################################
############################### Univariate Analysis ###############################
###################################################################################

par(mfrow=c(1,1))
plot.ts(full_stock_ts, main = "Microsoft closing monthly stock price (full)")
plot.ts(stock_ts, main = "Microsoft closing monthly stock price (short)")
plot.ts(infl_expect_ts, main = "10-YR Breakeven inflation rate")

############################## Stationarity & Trend ###############################

# Unit root test to check stationarity of the infl_expect
library(CADFtest)
max.lag<-round(sqrt(length(infl_expect_ts)))
CADFtest(infl_expect_ts, type= "drift", criterion= "BIC", max.lag.y=max.lag)
#p-value < 0.05 RH0 -> stationary (no unit root)

############ Trend ################################################################

# Trend for stock_ts
log_stock_ts <- log(stock_ts)
par(mfrow=c(1,1))
plot.ts(log_stock_ts)
TREND <- 1:203
fit <- lm(log_stock_ts ~ TREND)
summary(fit) #RH0 the trend is significant.
ts.plot(fit$residuals)
#Unit root test for the residuals
max2.lag<-round(sqrt(length(log_stock_ts)))
CADFtest(fit$residuals, type= "drift", criterion= "BIC", max.lag.y=max2.lag)
#Fail to reject H0 -> stochastic trend

# Trend for the full stock_ts is clearly stochastic
par(mfrow=c(1,1))
plot.ts(log(full_stock_ts))




############ Stationarity #############################################################

##### For stock_ts (Short time series) ###########################

# Going into log-differences for the stock_ts to make it stationary
dlog_stock_ts <- diff(log(stock_ts))
par(mfrow=c(1,1))
plot.ts(dlog_stock_ts, main = "Log-diff short MSFT stock price") #looks stationary



max5.lag<-round(sqrt(length(dlog_stock_ts)))

par(mfrow=c(1,1))
plot.ts(log_stock_ts)

par(mfrow=c(2,1))
acf(log_stock_ts)
pacf(log_stock_ts)

par(mfrow=c(2,1))
acf(stock_ts)
pacf(stock_ts)

par(mfrow=c(1,2))
acf(dlog_stock_ts)
pacf(dlog_stock_ts)

#White noise
Box.test(dlog_stock_ts, lag = max5.lag, type = "Ljung-Box") #Fail to reject H0: white noise
Box.test(infl_expect_ts, lag = max.lag, type = "Ljung-Box")



######################### Forceast models (short time series) ###############################


#1 log-diff, ARIMA(1,0,0)
mymodel1 <- arima(dlog_stock_ts, order=c(1,0,0))
mymodel1

par(mfrow=c(1,1))
plot.ts(mymodel1$res)
par(mfrow=c(1,2))
acf(mymodel1$res, main="Correlgram of ARIMA (1,0,0) residuals")
acf(mymodel1$res^2, main="Correlgram of ARIMA (1,0,0) residuals^2")
Box.test(mymodel1$res, lag = max2.lag, type = "Ljung-Box") #Fail to reject H0: white noise



#2 log-diff, ARIMA(0,0,1)
mymodel2 <- arima(dlog_stock_ts, order=c(0,0,1))
mymodel2

par(mfrow=c(1,1))
plot.ts(mymodel2$res)
par(mfrow=c(1,2))
acf(mymodel2$res, main="Correlgram of ARIMA (0,0,1) residuals")
acf(mymodel2$res^2, main="Correlgram of ARIMA (0,0,1) residuals^2")
Box.test(mymodel2$res, lag = max2.lag, type = "Ljung-Box") #Fail to reject H0: white noise


#3 log-diff, ARIMA(1,0,1)
mymodel3 <- arima(dlog_stock_ts, order=c(1,0,1))
mymodel3

par(mfrow=c(1,1))
plot.ts(mymodel3$res)
par(mfrow=c(1,2))
acf(mymodel3$res, plot = T, main="Correlgram of ARIMA (1,0,1) residuals")
acf(mymodel3$res^2, plot = T, main="Correlgram of ARIMA (1,0,1) residuals^2")
Box.test(mymodel3$res, lag = max2.lag, type = "Ljung-Box") #Fail to reject H0: white noise


#4 log-diff, ARIMA(2,0,1)
mymodel4 <- arima(dlog_stock_ts, order=c(2,0,1))
mymodel4

par(mfrow=c(1,1))
plot.ts(mymodel4$res)
par(mfrow=c(1,2))
acf(mymodel4$res, main="Correlgram of ARIMA (2,1,1) residuals")
acf(mymodel4$res^2, main="Correlgram of ARIMA (2,1,1) residuals^2")
par(mfrow=c(1,1))
pacf(mymodel4$res, main="Partial correlgram of ARIMA (2,1,1) residuals")
Box.test(mymodel4$res, lag = max2.lag, type = "Ljung-Box") #Fail to reject H0: white noise




#5 Auto ARIMA #########################
library(forecast)

proposed1 <- auto.arima(log_stock_ts) #ARIMA(0,1,0) with drift
proposed1
mymodel9 <- Arima(log_stock_ts,order = c(0,1,0), include.drift = TRUE)
mymodel9

par(mfrow=c(1,2))
acf(mymodel9$res, main="Correlgram for ARIMA(0,1,0) with drift res")
acf(mymodel9$res^2, main="Correlgram for ARIMA(0,1,0) with drift res^2")
par(mfrow=c(1,1))
pacf(mymodel9$res, main="Partial correlgram for ARIMA(0,1,0) with drift residuals")
Box.test(mymodel9$res, lag = max2.lag, type = "Ljung-Box") #Fail to reject H0: white noise



##################################################################


# Test the full series in log-differences for white noise


log_full_stock_ts <- log(full_stock_ts)
max3.lag<-round(sqrt(length(log_full_stock_ts)))

dlog_full_stock_ts <- diff(log_full_stock_ts)
max4.lag<-round(sqrt(length(dlog_full_stock_ts)))

Box.test(dlog_full_stock_ts, lag = max4.lag, type = "Ljung-Box") #white noise



par(mfrow=c(1,1))
plot.ts(log_full_stock_ts)
plot.ts(dlog_full_stock_ts, main = "Log-diff full MSFT stock price") #Stationary
par(mfrow=c(1,2))
acf(dlog_full_stock_ts,main="Correlgram for log-differences of the full stock series")
pacf(dlog_full_stock_ts,main="Partial correlgram for log-differences of the full stock series")

#1 ARIMA (1,0,1)
dmymodelf1 <- arima(dlog_full_stock_ts, order=c(1,0,1))
dmymodelf1

par(mfrow=c(1,1))
pacf(dmymodelf1$res, main="Partial correlgram of ARIMA (1,0,1) residuals")
par(mfrow=c(1,2))
acf(dmymodelf1$res, main="Correlgram of ARIMA (1,0,1) residuals")
acf(dmymodelf1$res^2, main="Correlgram of ARIMA (1,0,1) residuals^2")
Box.test(dmymodelf1$res, lag = max4.lag, type = "Ljung-Box") #Fail to reject H0: white noise

#1 ARIMA (2,0,1)
dmymodelf2 <- arima(dlog_full_stock_ts, order=c(2,0,1))
dmymodelf2

par(mfrow=c(1,1))
pacf(dmymodelf2$res, main="Partial correlgram of ARIMA (2,0,1) residuals")
par(mfrow=c(1,2))
acf(dmymodelf2$res, main="Correlgram of ARIMA (2,0,1) residuals")
acf(dmymodelf2$res^2, main="Correlgram of ARIMA (2,0,1) residuals^2")
Box.test(dmymodelf2$res, lag = max4.lag, type = "Ljung-Box") #Fail to reject H0: white noise

#2 ARIMA (1,1,2)
mymodelf2 <- arima(log_full_stock_ts, order=c(1,1,2))
mymodelf2

par(mfrow=c(1,1))
pacf(mymodelf2$res, main="Partial correlgram of ARIMA (1,1,2) residuals")
par(mfrow=c(1,2))
acf(mymodelf2$res, main="Correlgram of ARIMA (1,1,2) residuals")
acf(mymodelf2$res^2, main="Correlgram of ARIMA (1,1,2) residuals^2")
Box.test(mymodelf2$res, lag = max3.lag, type = "Ljung-Box") #Fail to reject H0: white noise

#3 ARIMA (3,1,0)
mymodelf3 <- arima(log_full_stock_ts, order=c(3,1,0))
mymodelf3

par(mfrow=c(1,1))
pacf(mymodelf3$res, main="Partial correlgram of ARIMA (3,1,0) residuals")
par(mfrow=c(1,2))
acf(mymodelf3$res, main="Correlgram of ARIMA (3,1,0) residuals")
acf(mymodelf3$res^2, main="Correlgram of ARIMA (3,1,0) residuals^2")
Box.test(mymodelf3$res, lag = max3.lag, type = "Ljung-Box") #Fail to reject H0: white noise

#4 ARIMA (0,1,3)
mymodelf4 <- arima(log_full_stock_ts, order=c(0,1,3))
mymodelf4

par(mfrow=c(1,1))
pacf(mymodelf4$res, main="Partial correlgram of ARIMA (0,1,3) residuals")
par(mfrow=c(1,2))
acf(mymodelf4$res, main="Correlgram of ARIMA (0,1,3) residuals")
acf(mymodelf4$res^2, main="Correlgram of ARIMA (0,1,3) residuals^2")
Box.test(mymodelf4$res, lag = max3.lag, type = "Ljung-Box") #Fail to reject H0: white noise



###########################################################

# Possible monthly effect
#Linear regression
lr_stock_ts = window(full_stock_ts, start = c(1986,3), end = c(2019,3))
dlog_lr_stock_ts = diff(log(lr_stock_ts))

par(mfrow=c(1,1))
monthplot(dlog_lr_stock_ts)

trend2 <- 1:396
M1<-rep(c(1,0,0,0,0,0,0,0,0,0,0,0),33)
M2<-rep(c(0,1,0,0,0,0,0,0,0,0,0,0),33)
M3<-rep(c(0,0,1,0,0,0,0,0,0,0,0,0),33)
M4<-rep(c(0,0,0,1,0,0,0,0,0,0,0,0),33)
M5<-rep(c(0,0,0,0,1,0,0,0,0,0,0,0),33)
M6<-rep(c(0,0,0,0,0,1,0,0,0,0,0,0),33)
M7<-rep(c(0,0,0,0,0,0,1,0,0,0,0,0),33)
M8<-rep(c(0,0,0,0,0,0,0,1,0,0,0,0),33)
M9<-rep(c(0,0,0,0,0,0,0,0,1,0,0,0),33)
M10<-rep(c(0,0,0,0,0,0,0,0,0,1,0,0),33)
M11<-rep(c(0,0,0,0,0,0,0,0,0,0,1,0),33)
M12<-rep(c(0,0,0,0,0,0,0,0,0,0,0,1),33)
fit2<-lm(dlog_lr_stock_ts~trend2 + M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11+M12)
summary(fit2)
# No significant monthly effect, model not significant

#####################################################################

#Add GARCH to handle Heteroscedasticity
library(fGarch)
#Data needs to be stationary (use serie in log-diff)

#ARMA(1,1)+GARCH(1,1)
fit_garch0<-garchFit(~arma(1,1)+garch(1,1),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch0
summary(fit_garch0)
par(mfrow=c(1,1))
#plot(fit_garch0)


#Out of sample expanding Window MAE calculation
y<-dlog_full_stock_ts
S=round(0.75*length(y))
h=1
error1.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(1,1)+garch(1,1),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
MAE1<-mean(abs(error1.h))
MAE1


#ARMA(1,1)+GARCH(2,2)
fit_garch5<-garchFit(~arma(1,1)+garch(2,2),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch5
summary(fit_garch5)
par(mfrow=c(1,1))
#plot(fit_garch5)

#Out of sample expanding Window MAE calculation
error2.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(1,1)+garch(2,2),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
MAE2<-mean(abs(error2.h))
MAE2


#ARMA(2,1)+GARCH(1,1)
fit_garch<-garchFit(~arma(2,1)+garch(1,1),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch
summary(fit_garch)
par(mfrow=c(1,1))
#plot(fit_garch)

#Out of sample expanding Window MAE calculation
error3.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(2,1)+garch(1,1),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error3.h<-c(error3.h,y[i+h]-predict.h)
}
MAE3<-mean(abs(error3.h))
MAE3


#ARMA(1,2)+GARCH(1,1)
fit_garch2<-garchFit(~arma(1,2)+garch(1,1),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch2
summary(fit_garch2)
par(mfrow=c(1,1))
#plot(fit_garch2)

#Out of sample expanding Window MAE calculation
error4.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(1,2)+garch(1,1),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error4.h<-c(error4.h,y[i+h]-predict.h)
}
MAE4<-mean(abs(error4.h))
MAE4


#ARMA(3,0)+GARCH(1,1)
fit_garch3<-garchFit(~arma(3,0)+garch(1,1),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch3
summary(fit_garch3)
par(mfrow=c(1,1))
#plot(fit_garch3)

#Out of sample expanding Window MAE calculation
error5.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(3,0)+garch(1,1),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error5.h<-c(error5.h,y[i+h]-predict.h)
}
MAE5<-mean(abs(error5.h))
MAE5


#ARMA(0,3)+GARCH(1,1)
fit_garch4<-garchFit(~arma(0,3)+garch(1,1),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch4
summary(fit_garch4)
par(mfrow=c(1,1))
#plot(fit_garch4)

#Out of sample expanding Window MAE calculation
error6.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(0,3)+garch(1,1),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error6.h<-c(error6.h,y[i+h]-predict.h)
}
MAE6<-mean(abs(error6.h))
MAE6


#ARMA(2,2)+GARCH(1,1)
fit_garch6<-garchFit(~arma(2,2)+garch(1,1),dist="QMLE",include.mean=T, data=dlog_full_stock_ts) 
fit_garch6
summary(fit_garch6)
par(mfrow=c(1,1))
#plot(fit_garch6)

#Out of sample expanding Window MAE calculation
error7.h<-c()
for (i in S:(length(y)-h))
{
  mymodel.sub<-garchFit(~arma(2,2)+garch(1,1),dist="QMLE",include.mean=T, data=y[1:i])
  predict.h<-predict(mymodel.sub,n.ahead=h)$meanForecast[h]
  error7.h<-c(error7.h,y[i+h]-predict.h)
}
MAE7<-mean(abs(error7.h))
MAE7

#####################################################
### Auto.rima For the full model

library(forecast)

proposed2 <- auto.arima(log_full_stock_ts) #SARMA (1,2,1)(2,0,0)
proposed2
mymodelf9 <- arima(log_full_stock_ts,order = c(1,2,1), seasonal = c(2,0,0))
mymodelf9

par(mfrow=c(1,2))
acf(mymodelf9$res, main="Correlgram for SARMA (1,2,1)(2,0,0) residuals")
acf(mymodelf9$res^2, main="Correlgram for SARMA (1,2,1)(2,0,0) residuals^2")
par(mfrow=c(1,1))
pacf(mymodelf9$res, main="Partial correlgram for SARMA (1,2,1)(2,0,0) residuals")
Box.test(mymodelf9$res, lag = max3.lag, type = "Ljung-Box") #Fail to reject H0: white noise
#Model rejected due to Heteroscedasticity

proposed3 <- auto.arima(dlog_full_stock_ts) #SARMA (2,1,3)(2,0,0)
proposed3
mymodelf10 <- arima(dlog_full_stock_ts,order = c(2,1,3), seasonal = c(2,0,0))
mymodelf10

par(mfrow=c(1,2))
acf(mymodelf10$res, main="Correlgram for SARMA (1,2,1)(2,0,0) residuals")
acf(mymodelf10$res^2, main="Correlgram for SARMA (1,2,1)(2,0,0) residuals^2")
par(mfrow=c(1,1))
pacf(mymodelf10$res, main="Partial correlgram for SARMA (1,2,1)(2,0,0) residuals")
Box.test(mymodelf10$res, lag = max4.lag, type = "Ljung-Box") #Fail to reject H0: white noise
#Model rejected due to Heteroscedasticity

#########################################################


#Diebold-Mariano test between ARMA (1,1) + GARCH (1,1) and ARMA (1,2) + GARCH (1,1) 
#to determine significant difference

library(forecast)
dm.test(error1.h,error4.h,h=h,power=1)# Modified Diebold Mariano test (1997) alrready has a correction for autocorrelation
#Fail to reject H0: there is no significant difference

# ARMA (1,1) + GARCH (1,1) model is chosen since it is simpler.

# Forecast
library(fGarch)

myforecast<-predict(fit_garch0,n.ahead=12, plot = T)
myforecast


#####################################################################################
############################### Multivariate Analysis ###############################
#####################################################################################

# Note that the 10-Year Breakeven Inflation rate, is already in differences (rate) and stationary
# There is a limitation for some calculations like Engle-Granger test since the 10-Year Breakeven Inflation Rate is I(0).
# The corresponding series I(1) will be calculated based on US December-2002 CPI of 180.9 (1982-84=100).

CPI <- c()
CPI[1] <- 180.9

for (i in 2:203)
{
  CPI[i]<- CPI[i-1]*(1 + infl_expect_ts[i-1]/100)
}

CPI_ts <- ts(CPI, start = c(2002,12), frequency = 12)

logCPI_ts <- log(CPI_ts)*100
plot.ts(logCPI_ts)
plot.ts(diff(logCPI_ts))


# ADL(3)
lag <- 3
n <- length(dlog_stock_ts)
dlog_CPI_ts <- diff(logCPI_ts)
dlogStock.0 <- dlog_stock_ts[(lag+1):n]
dlogStock.1 <- dlog_stock_ts[lag:(n-1)]
Infl.1 <- dlog_CPI_ts[lag:(n-1)]
Infl.2 <- dlog_CPI_ts[(lag-1):(n-2)]
Infl.3 <- dlog_CPI_ts[(lag-2):(n-3)]
fit_dlm <- lm(dlogStock.0 ~ dlogStock.1+Infl.1+Infl.2+Infl.3)
summary(fit_dlm) #RH0: model is significant (not very significant)
par(mfrow=c(1,1))
acf(fit_dlm$residuals) #Border line significance at lag 4 and 15
Box.test(fit_dlm$residuals, lag = max5.lag, type = "Ljung-Box") #Fail to reject H0: white noise

# ADL(4)
lag2 <- 4
dlogStock2.0 <- dlog_stock_ts[(lag2+1):n]
dlogStock2.1 <- dlog_stock_ts[lag2:(n-1)]
Infl2.1 <- dlog_CPI_ts[lag2:(n-1)]
Infl2.2 <- dlog_CPI_ts[(lag2-1):(n-2)]
Infl2.3 <- dlog_CPI_ts[(lag2-2):(n-3)]
Infl2.4 <- dlog_CPI_ts[(lag2-3):(n-4)]
fit_dlm2 <- lm(dlogStock2.0 ~ dlogStock2.1+Infl2.1+Infl2.2+Infl2.3+Infl2.4)
summary(fit_dlm2) #RH0: model is significant (not very significant)
par(mfrow=c(1,1))
acf(fit_dlm2$residuals) #Border line significance at lag 4 and 15
Box.test(fit_dlm2$residuals, lag = max5.lag, type = "Ljung-Box") #Fail to reject H0: white noise

fit_dlm

# Granger Causality 

fit_dlm_nox <- lm(dlogStock.0 ~ dlogStock.1)
anova(fit_dlm,fit_dlm_nox)
# RH0 of no Granger Causality: There is a Granger Causality. 
# 10-Year Breakeven Inflation Rate (lag 3) has an incremental power in predictiong dlog MSFT Stock


# Engle-Granger test for no cointegration
library(CADFtest)

fit_ci<-lm(log_stock_ts ~ logCPI_ts)
res_fit_ci<-fit_ci$residuals
CADFtest(res_fit_ci,type="drift",criterion="BIC",max.lag.y=max2.lag)
# Fail to reject H0. T-statstic (-0.78165) is greater than -3.41 => H0: No evidence for cointegration (No cointegration)


# Var model with the stationary series (No ECM nor VECM)

library(vars)
dlogdata<-data.frame(dlog_stock_ts,dlog_CPI_ts)

VARselect(dlogdata,lag.max=10,type="const") # Order selected 2
fit_varautom<-VAR(dlogdata,type="const",p=2)
summary(fit_varautom)
varautom_residuals<-resid(fit_varautom)

par(mfrow=c(2,2))
acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
ccf(varautom_residuals[,1],varautom_residuals[,2])

# Impulse response
irf_var<-irf(fit_varautom,ortho=F,boot=T)
plot(irf_var)

# Var model is not a good approximation for dlog_stock_ts

# There is no significance for dlog_stock_ts = dlog_stock_ts.l1 + dlog_CPI_ts.l1 + dlog_stock_ts.l2 + dlog_CPI_ts.l2 + const 
# Some estimators are significant, but not the whole model.
# However dlog_CPI_ts = dlog_stock_ts.l1 + dlog_CPI_ts.l1 + dlog_stock_ts.l2 + dlog_CPI_ts.l2 + const  is significant.













