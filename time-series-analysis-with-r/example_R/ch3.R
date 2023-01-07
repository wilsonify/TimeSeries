library(astsa)
library(TSA)


#Section 3.3: ACF and PACF 

#AUTOREGRESSIVE(p)

# Sample Paths of AR(1)
# in the expressions below, ~ is a space and == is equal
ar11=arima.sim(list(order=c(1,0,0), ar=.9), n=100)
ar12=arima.sim(list(order=c(1,0,0), ar=.4), n=100)
ar13=arima.sim(list(order=c(1,0,0), ar=-.9), n=100)
ar14=arima.sim(list(order=c(1,0,0), ar=-.5), n=100)

par(mfrow=c(2,2))     
ts.plot(ar11, ylab="x", main=(expression(AR(1)~~~phi==+.9))) 
ts.plot(ar12, ylab="x", main=(expression(AR(1)~~~phi==+.4))) 
ts.plot(ar13, ylab="x", main=(expression(AR(1)~~~phi==-.9))) 
ts.plot(ar14, ylab="x", main=(expression(AR(1)~~~phi==-.5))) 

 



#True ACFs of AR(1)
h=2:20
par(mfrow=c(2,2))                         
plot(h,(0.9)^h, type="h", main=(expression(AR(1)~~~phi==+.9)) )
points(h,(0.9)^h)
abline(h=0)
plot(h,(0.4)^h, type="h", main=(expression(AR(1)~~~phi==+.4)) )
points(h,(0.4)^h)
abline(h=0)
plot(h,(-0.9)^h, type="h", main=(expression(AR(1)~~~phi==-.9)) )
points(h,(-0.9)^h)
abline(h=0)
plot(h,(-0.5)^h, type="h", main=(expression(AR(1)~~~phi==-.5)) )
points(h,(-0.5)^h)
abline(h=0)


#Sample ACF's
par(mfrow=c(2,2)) 
acf(ar11,20,main=(expression(Simu~AR(1)~~~phi==+.9)))
acf(ar12,20,main=(expression(Simu~AR(1)~~~phi==+.4)))
acf(ar13,20,main=(expression(Simu~AR(1)~~~phi==-.9)))
acf(ar14,20,main=(expression(Simu~AR(1)~~~phi==-.5)))



#AR(2) Simulation
par(mfrow=c(2,2))
ar21 <- arima.sim(list(order = c(2,0,0), ar = c(0.5, 0.25)), n = 100)
ar22 <- arima.sim(list(order = c(2,0,0), ar = c(1,-0.25)), n = 100)
ar23 <- arima.sim(list(order = c(2,0,0), ar = c(1.5,-0.75)), n = 100)
ar24 <- arima.sim(list(order = c(2,0,0), ar = c(1,-0.6)), n = 100)
plot(ar21,ylab=expression(x[t]),xlab="Time",type="o")
plot(ar22,ylab=expression(x[t]),xlab="Time",type="o")
plot(ar23,ylab=expression(x[t]),xlab="Time",type="o")
plot(ar24,ylab=expression(x[t]),xlab="Time",type="o")

#Root calculation  determine stationarity of AR(2)
#Be careful with the signs as you formulate the characteristic equation       

#i)
#polynomial formulation in R:  1+(-0.5)*B+(-0.25)*B^2
z = c(1,-0.5, -.25)   
Mod( polyroot(z)   ) #modulus of a complex number

(a = Mod(polyroot(z)[1])) #  print one root  
 

#ii)
z = c(1,-1,.25) #polynomial: 1+(-1)*B +(0.25)*B^2
a=polyroot(z)
a
Mod(a) #outside the unit circle BUT has only one root = 2  
 

#iii)
z = c(1,-1.5,.75)  #polynomial: 1+(-1.5)*B +(0.75)*B^2
    
polyroot(z)
(a = polyroot(z)[1])   
arg = Arg(a)/(2*pi)  # arg in cycles/unit time   (frequency)
1/arg                # = 12,  the period  
Mod(polyroot(z))
  

#iv)
z = c(1,-1,.6)    #polynomial: 1+(-1)*B +(0.6)*B^2
polyroot(z) 
Mod(polyroot(z))
(a = polyroot(z)[1])  
arg = Arg(a)/(2*pi)  # arg in cycles/unit time   (frequency)
1/arg                # = 7,  the period 
 

#########################
#True ACF of AR(2).  Uses ARMAacf function of TSA
# ARMAacf function includes the k=0 lag
# Use y = y[2:21] to remove k=0 lag from ARMAacf output
library(TSA)
par(mfrow=c(2,2))
y = ARMAacf(ar = c(0.5,0.25), lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "h",
      ylab = "Autocorrelation", main = "Population ACF")
abline(h = 0)
y = ARMAacf(ar = c(1,-0.25), lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "h",
      ylab = "Autocorrelation", main = "Population ACF")
abline(h = 0)
y = ARMAacf(ar = c(1.5,-0.75), lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "h",
      ylab = "Autocorrelation", main = "Population ACF")
abline(h = 0)
y = ARMAacf(ar = c(1,-0.6), lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "h",
      ylab = "Autocorrelation", main = "Population ACF")
abline(h = 0)
 



# AR(2) sample ACFs
par(mfrow=c(2,2))
acf(ar21,main="Sample ACF")
acf(ar22,main="Sample ACF")
acf(ar23,main="Sample ACF")
acf(ar24,main="Sample ACF")

 
#AR(p) Applications

#Monthly deaths from bronchitis, emphysema and asthma in the UK
#1974–1979, both sexes (ldeaths), males (mdeaths) and females (fdeaths).
 
library(datasets)
par(mfrow=c(3,1))
plot(ldeaths)
plot(mdeaths)
plot(fdeaths)

par(mfrow=c(3,1))
acf(ldeaths, lag.max=70)
acf(mdeaths, lag.max=70)
acf(fdeaths,lag.max=70)
 


###############################################################
#MOVING AVERAGE (q)
 

#MA(1)
ma11=arima.sim(list(order=c(0,0,1), ma=.9), n=100)                              
ma12=arima.sim(list(order=c(0,0,1), ma= -.9), n=100)  
                            
par(mfrow=c(2,1)) 
ts.plot(ma11, ylab="x", main=(expression(MA(1)~~~theta==+.9 )))    
ts.plot(ma12, ylab="x", main=(expression(MA(1)~~~theta==-.9)))   

#Sample ACF's
par(mfrow=c(2,1)) 
acf(ma11,20,main=(expression(Simu~MA(1)~~~theta==+.9)))
acf(ma12,20,main=(expression(Simu~MA(1)~~~theta==-.9)))


#Non-uniqueness of MA(1)
thet=seq(-1,1,length.out=100)
plot(thet,thet/(1+ thet^2) )
points(thet,(1/thet)/(1+ (1/thet)^2), col="red" )
lines(thet,thet/(1+ thet^2))
lines(thet,(1/thet)/(1+ (1/thet)^2), col="red" )

#Scatter Plots for simulated MA(1)  data
library(astsa)

#lag-1 scatter plot
ma11b=ma11[2:length(ma11)] # lag(as.matrix(ma11), -1, na.pad = TRUE), lag is not working-not sure
ma12b=ma12[2:length(ma12)] 

par(mfrow=c(1,2))
plot( ma11b,   as.matrix(ma11)[1:(length(ma11)-1)],type="p",xlab=expression(x[t-1]),ylab=expression(x[t]),main="Lag 1 scatterplot, theta==+.9" )
plot( ma12b,   as.matrix(ma12)[1:(length(ma12)-1)],type="p",xlab=expression(x[t-1]),ylab=expression(x[t]),main="Lag 1 scatterplot, theta==-.9")
 
 
#lag-2 scatter plot
ma11c=ma11[3:length(ma11)] # lag(as.matrix(ma11), -1, na.pad = TRUE), lag is not working-not sure
ma12c=ma12[3:length(ma12)] 

par(mfrow=c(1,2))
plot( ma11c,   as.matrix(ma11)[1:(length(ma11)-2)],type="p",xlab=expression(x[t-1]),ylab=expression(x[t]),main="Lag 2 scatterplot, theta==+.9" )
plot( ma12c,   as.matrix(ma12)[1:(length(ma12)-2)],type="p",xlab=expression(x[t-1]),ylab=expression(x[t]),main="Lag 2 scatterplot, theta==-.9")
 
 

##############
#MA(2) 
# An MA(2) series with MA coefficients equal to 1 and 0.9 and 
# of length n=100 can be simulated by the following command

ma21=arima.sim(model=list(ma= c(1, 0.9)),n=100) 
#ma21=arima.sim(list(order=c(0,0,2), ma= c(1, 0.9) ), n=100)

library(TSA) #for zlag to work
par(mfrow=c(2,2))
plot(ma21,ylab=expression(x[t]),xlab="Time",type="o",main="MA(2) simulation")
acf(ma21,main="Sample ACF")$acf
plot(y=ma21, x=zlag(ma21,1),ylab=expression(x[t]),xlab=expression(x[t-1]),type='p',main="Lag 1 scatterplot")
plot(y=ma21,x=zlag(ma21,2),ylab=expression(x[t]),xlab=expression(x[t-2]),type='p',main="Lag 2 scatterplot")

#rho1
(1+1*0.9)/(1+1^1+0.9^2)

#rho2
0.9/(1+1^1+0.9^2)


#Invertibility of MA(2)
z = c(1,1,.9)
Mod(polyroot(z))    #>1
 

#Invertibility of MA(2): Another example
z = c(1,-2,2)
Mod(polyroot(z))    #<1
 

#US CONSUMPTION APPLICATION
library(datasets)
library(fpp)
data(usconsumption)
dim(usconsumption)
par(mfrow=c(2,1))
plot(usconsumption[,1])
acf(usconsumption[,1])

 
 
#######################################
#######################################
#ARMA(p,q)
#######################################
#######################################

#See slides
#ARMA(p,q) simulation
par( mfrow=c(2,2) )
#ARMA(1,1)
arma11.sim <- arima.sim(list(order = c(1,0,1), ar = 0.6, ma = 0.8), n = 200)
ts.plot(arma11.sim) 

#ARMA(2,1)
arma21.sim <- arima.sim(list(order = c(2,0,1), ar = c(-0.3, 0.6), ma = 0.8), n = 200)
ts.plot(arma21.sim) 


#ARMA(2,2)
arma22.sim <- arima.sim(list(order = c(2,0,2), ar = c(0.5,-0.5), ma = c(0.8,-0.2)), n = 200)
ts.plot(arma22.sim)  

#ARMA(3,3)
arma33.sim <- arima.sim(list(order = c(3,0,3), ar = c(0.8,0.8,-0.9), ma = c(-0.9,0.8,-0.2)), n = 200)
ts.plot(arma33.sim) 


#########################
#TRUE ARMA ACF 
#########################

# True ARMA(1,1) autocorrelation functions 
# Uses ARMAacf function
# ARMAacf function includes the k=0 lag
# Use y = y[2:21] to remove k=0 lag from ARMAacf output
# Page 112
par(mfrow=c(2,2))
y = ARMAacf(ar = 0.8, ma = 0.2, lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Autocorrelation", main = "True ACF")
abline(h = 0)
y = ARMAacf(ar = -0.8, ma = 0.2, lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Autocorrelation", main = "True ACF")
abline(h = 0)
y = ARMAacf(ar = 0.6, ma = -0.3, lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Autocorrelation", main = "True ACF")
abline(h = 0)
y = ARMAacf(ar = -0.6, ma = -0.3, lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Autocorrelation", main = "True ACF")
abline(h = 0)

##############
#ARMA(p,q) Application/s 
############
library(tseries)
data(nino)
s <- nino3.4
par(mfrow=c(2,1))
plot(s)
acf(s, lag.max=40)  
 

#Monthly deaths from bronchitis, emphysema and asthma in the UK
#1974–1979, both sexes (ldeaths), males (mdeaths) and females (fdeaths).

acf(ldeaths, lag.max=70)
acf(mdeaths, lag.max=70)
acf(fdeaths,lag.max=70)


#############################
#############################
#PACF
#############################
#############################

library(astsa)

#Example 3.16: PACF of AR(1) and AR(2)
ar1.acf = ARMAacf(ar=0.5, ma=0, 24)   
ar1.pacf = ARMAacf(ar=0.5, ma=0, 24, pacf=TRUE) #excludes h=0
ar2.acf = ARMAacf(ar=c(1.5,-.75), ma=0, 24) 
ar2.pacf = ARMAacf(ar=c(1.5,-.75), ma=0, 24, pacf=TRUE)#excludes h=0

par(mfrow=c(2,2))
plot(ar1.acf, type="h", xlab="lag")
abline(h=0)
plot(ar1.pacf, type="h", xlab="lag")
abline(h=0)
plot(ar2.acf, type="h", xlab="lag")
abline(h=0)
plot(ar2.pacf, type="h", xlab="lag")
abline(h=0)

# AR(1) and AR(2) simulated ACFs and PACFs
ar1.sim <- arima.sim(list(order = c(1,0,0), ar = c(0.5)), n = 150)
ar2.sim <- arima.sim(list(order = c(2,0,0), ar = c(1.5,-0.75)), n = 150)
plot(ar1.sim,ylab=expression(x[t]),main="AR(1) process",xlab="Time",type="o")
plot(ar2.sim,ylab=expression(x[t]),main="AR(2) process",xlab="Time",type="o")
acf2(ar1.sim) #use acf and pacf in TSA
acf2(ar2.sim) #acf2 from astsa package



# MA(1) and MA(2) TRUE ACF/PACF
par(mfrow=c(2,2))
y = ARMAacf(ma = -0.8, lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Autocorrelation", main = "True ACF")
abline(h = 0)
y = ARMAacf(ma = -0.8, lag.max = 20, pacf=TRUE)
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Partial autocorrelation", main = "True PACF")
abline(h = 0)
y = ARMAacf(ma = c(0.6,-0.3), lag.max = 20)
y = y[2:21]
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Autocorrelation", main = "True ACF")
abline(h = 0)
y = ARMAacf(ma = c(0.6,-0.3), lag.max = 20, pacf=TRUE)
y
plot(y, x = 1:20, type = "h", ylim = c(-1,1), xlab = "k",
      ylab = "Partial autocorrelation", main = "True PACF")
abline(h = 0)


library(astsa)
#El Nino's Recruitment data
acf2(rec, 24)

#Luteinizing Hormone
acf2(lh,24)

#monthly deaths from bronchitis, emphysema and asthma in the UK, 1974–1979, 
#both sexes (ldeaths)
acf2(ldeaths) 
 

#####################################
#####################################
#Extended ACFs
#####################################
#####################################


library(astsa)
library(TSA)
#ARMA(p,q) simulations, ACF,PACF
#Sample EACF calculation, pp 116, CaC
#ARMA(1,1)
arma11.sim <- arima.sim(list(order = c(1,0,1), ar = 0.6, ma = 0.8), n = 200)
acf2(arma11.sim) #ACF and PACF plots
eacf(arma11.sim)$eacf #estimates
eacf(arma11.sim)$symbol #limits

#ARMA(2,2)
arma22.sim <- arima.sim(list(order = c(2,0,2), ar = c(0.5,-0.5), ma = c(0.8,-0.2)), n = 200)
acf2(arma22.sim) #ACF and PACF plots
eacf(arma22.sim)$eacf #estimates
eacf(arma22.sim)$symbol #limits


#ARMA(3,3)
arma33.sim <- arima.sim(list(order = c(3,0,3), ar = c(0.8,0.8,-0.9), ma = c(-0.9,0.8,-0.2)), n = 200)
acf2(arma33.sim) #ACF and PACF plots
eacf(arma33.sim)$eacf #estimates
eacf(arma33.sim)$symbol #limits
 

#AUTO.ARIMA
library(forecast)
#install.packages("forecast")
#install.packages("xts")
auto.arima(arma33.sim)
 

##############################
#More APPLICATIONS

#Monthly Yields On Treasury Securities
This data set contains monthly 1 year, 3 year, 
#5 year, and 10 year yields on treasury securities at constant, fixed maturity.
#four univariate time series tcm1y, tcm3y, tcm5y, and tcm10y and the joint series tcm.
#from 1953-1999

#install.packages("tseries")
#install.packages("TTR")
library(tseries)

data(tcm)  
par(mfrow=c(2,1))
plot(tcm10y)
acf(tcm10y)

library(astsa)
dat <- diff(tcm10y)
acf2(dat) 
eacf(dat)

library(forecast)
auto.arima(dat)

#########################
#More Exercises   
#Investigate the ACF and PACF of the El Nino's SOI  and Recruitmentdata.
  

########################################################
#3.5: ESTIMATION
########################################################
library(astsa)
#Example 3.27(text). AR(2): xt=1.5x(t-1)-.75x(t-2) +wt
set.seed(8675309)
ar2.sim=arima.sim(list(order = c(2,0,0), ar = c(1.5,-0.75)), n = 144)
acf2(ar2.sim)


###########
#estimation with goodness-of-fit diagnostics
sarima(ar2.sim, 2,0,0) 
sarima(ar2.sim, 2,0,0, no.constant=TRUE) 
 
#Example 3.29 (text): MA(1), xt= wt + theta w(t-1) , theta=0.9.
set.seed(2)
ma1 = arima.sim(list(order = c(0,0,1), ma = 0.9), n = 150)

###########
#estimation with goodness-of-fit diagnostics
sarima(ma1, 0,0,1) 
sarima(ma1, 0,0,1, no.constant=TRUE) 

#ARMA(1,1) , phi =0.83, theta= -0.43.
set.seed(2786)
arma11.sim = arima.sim(list(order = c(1,0,1), ar=0.83,  ma = -0.43), n = 150)
 
###########
#estimation with goodness-of-fit diagnostics
sarima(arma11.sim, 1,0,1) 

#############
#Recruitment data
#estimation with  goodness-of-fit diagnostics
acf2(rec)
sarima(rec, 2,0,0) 
 
#auto.arima
auto.arima(rec)
sarima(rec,1,1,0,0,0,2,12, no.constant=TRUE)


#Monthly lung disease deaths for females (fdeaths)
#auto.arima
acf2(fdeaths)
sarima(fdeaths, 4,0,0 )

auto.arima(fdeaths)
sarima(fdeaths,0,0,0,2,1,0,12, no.constant=TRUE)
 


 

###########################################
#Section 3.4: FORECAST                    # 
###########################################

#Global warming data: temp deviations (1880-2015)
library(astsa)
data(globtemp)
plot(globtemp)
str(globtemp)
time=time(globtemp)
time2=time^2
time(globtemp)
dat=cbind(globtemp,time,time2)
mod1=lm(globtemp~time+time2,dat)
mod1
#plot fit
plot(time,globtemp, type='l',xlim=c(1880,2015), ylim=c(-0.5, 1),frame.plot=FALSE)
lines( 2.844e+02   -2.992e-01*time +7.860e-05*time2, lwd=3, col="red")

#residuals
qqnorm(resid(mod1),frame.plot=FALSE)
qqline(resid(mod1))

library(astsa)
sarima.for(rec,24,2,0,0) #with +/- 1,2*SE (prediction error bounds)



#OR
library(astsa)
sarima.for(fdeaths,24,0,0,0,2,1,0,12) #with +/- 1,2*SE (prediction error bounds)


######################################
######################################
#Integrated Models for Non-Stationary Data
#or ARIMA
######################################
######################################

#ARI(1,1)
library(astsa)
par(mfrow=c(2,2))
ari11.sim <- arima.sim(list(order = c(1,1,0), ar = 0.5), n = 150)
plot(ari11.sim,ylab=expression(x[t]),xlab="Time",type="o")
acf(ari11.sim,main="Sample ACF: ARI(1,1)")
plot(diff(ari11.sim),ylab=expression(x[t]-x[t-1]),xlab="Time",type="o")
acf(diff(ari11.sim),main="Sample ACF: 1st differences")
 

#######################################
#GNP DATA
data(gnp)
 
par(mfrow=c(1,2))
plot(gnp,frame.plot=FALSE)
plot(log(gnp),frame.plot=FALSE)
acf(gnp, 50, frame.plot=FALSE)  

gnpgr = diff(log(gnp))  # growth rate or log-return is a common transformation for finance or stock data
plot(gnpgr)
acf2(gnpgr, 24)  

#fitting # MA(2) 
sarima(gnpgr, 0, 0, 2)      

 
#fitting  AR(1) 
sarima(gnpgr, 1, 0, 0)    
sarima(log(gnp), 1, 1, 0)   # we need the constant mu  
 

#OR
library(forecast)
auto.arima(log(gnp))

sarima(gnpgr, 1, 0, 0)   
sarima(log(gnp), 1, 1, 0) 

#FORECAST
library(astsa)
sarima.for(gnpgr,30,1,0,0 )
points(gnpgr,xlim=c(1970, 2010), col="blue")

#FORECASTING UNDIFFERENCED DATA
sarima.for(log(gnp),30,1,1,0 )
points(log(gnp),xlim=c(1970, 2010), col="blue")


  
######################################
#GLACIAL VALVE DATA
par(mfrow=c(1,2))
plot(varve,frame.plot=FALSE)
acf(varve, 50, frame.plot=FALSE)  

#log-transform
acf2(log(varve), 50, frame.plot=FALSE)  

#apply differencing
acf2(diff(log(varve)), 50, frame.plot=FALSE)  

####################################### 
#FITTING/ESTIMATION and DIAGNOSTICS

#BEWARE of "no.constant = TRUE" for ARIMA(p,d,q)
#if d=1, put no.constant = TRUE for ARIMA(p,d,q)
sarima(log(varve), 0, 1, 1) #constant is not significant
sarima(log(varve), 0, 1, 1, no.constant=TRUE)    # ARIMA(0,1,1) and "no.constant=TRUE" includes overall mean constant
sarima(log(varve), 1, 1, 1, no.constant=TRUE)   # ARIMA(1,1,1)

#OR
library(forecast)
auto.arima(log(varve))


#FORECASTING UNDIFFERENCED DATA
library(astsa) 
sarima.for( log(varve),30,1,1,1 )
points(log(varve),    col="blue")
 

#FORECASTING DIFFERENCED DATA
library(astsa) 
sarima.for( diff(log(varve)),30,1,0,1 )
points(diff(log(varve)),    col="blue")
 
##################################
#SARIMA
#################################
 
 
###############################
#SARIMA SIMULATION
#SARIMA(1,0,0) x (1,1,0)[12]
library(astsa)
library(CombMSC)
set.seed(375465)
sdat1<- sarima.Sim(n=10, period=12, model = list(order = c(1,0,0), ar=0.5),seasonal=list(order= c(1,1,0), ar = 0.8))
acf2(sdat1)
ggseasonplot(sdat1)
acf2(diff(sdat1,12))

#fit and goodness of fit
sarima(sdat1,1,0,0,1,1,0,12)

##############################
#Airline Passengers: Monthly total 
x     = AirPassengers
plot(x)

lx    = log(x) 
plot(lx)
 
dlx   = diff(lx)  

acf2(dlx)
 
#EDA to determine seasonality
library(forecast)
ggseasonplot(dlx)
 
ddlx  = diff(dlx, 12) #=(1-B^s)^D *x(t) = x(t)-x(t-s), D=1, s=12
acf2(ddlx,100) #determine p, q, P, Q 
eacf(ddlx)
#now stationary so we are ready to fit a model
plot.ts(cbind(x, lx, dlx, ddlx), main="")
 
sarima(lx, 1,1,1, 0,1,1, 12)   # model 1
sarima(lx, 0,1,1, 0,1,1, 12)   # model 2   
sarima(lx, 1,1,0, 0,1,1, 12)   # model 3

sarima.for(lx, 12, 0,1,1, 0,1,1,12)  # forecasts with +/- 1 and 2 prediction error bounds. 


#AUTO.ARIMA
auto.arima(lx)
sarima(lx, 0,1,1, 0,1,1, 12)   # same as model 2


##########################
#Monthly deaths from bronchitis, emphysema and asthma in the UK
#1974–1979, both sexes (ldeaths), males (mdeaths) and females (fdeaths).
  
acf2(ldeaths)  

ggseasonplot(ldeaths)

auto.arima(ldeaths) #with drift (constant)

sarima(ldeaths, 0,0,2, 2,1,0, 12) #auto.arima model
sarima(ldeaths, 2,0,0, 1,1,0,12)  #Box-jenkins
sarima(ldeaths, 4,0,0)  #Box-jenkins


#tweaking auto.arima via testing coeffcients
sarima(ldeaths, 0,0,1, 2,1,0, 12) 


#forecasting
sarima.for(ldeaths, 24, 0,0,1, 2,1,0, 12) 

##########################
#CHAPTER 4
##########################


#######
#SOI and REC MONTHLY data
par(mfrow=c(2,1))
acf(soi, lag.max=100)
acf(rec, lag.max=100)

#Smoothed Periodogram
par(mfrow=c(2,1))
soi.ave = mvspec(soi, kernel('daniell',4), log='no')
abline(v = c(.25,1,2,3), lty=2) 
# Repeat above commands with soi replaced by rec, for example:
rec.ave = mvspec(rec, kernel('daniell',4), log="no")
abline(v=c(.25,1,2,3), lty=2)
 
#RESULTS and INTEPRETATIONS
#The frequency axis  is labeled in multiples of 1/12.
#Most dominant period for SOI = 1/(1*(1/12)) = 12 months = 1 year
#Next dominant period for SOI = 1/( (1/4)*(1/12) )  = 48 months = 4 years 

#Most dominant period for REC = 1/( (1/4)*(1/12) )  = 48 months = 4 years 
#Next dominant period for REC =   1/(1*(1/12)) = 12 months = 1 year



 

