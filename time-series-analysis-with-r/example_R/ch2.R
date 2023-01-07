#Section 2.1
library(astsa)

#DETRENDING
#Chicken prices
ts.plot(chicken, ylab="cents per pound", col=4, lwd=2)
summary(fit <- lm(chicken~time(chicken))) # regress price on time
abline(fit)        # add the fitted regression line to the plot    
acf(resid(fit), 30) 
  

#Pollution, temperature and mortality for model (2.21) only
par(mfrow=c(3,1))
ts.plot(cmort, main="Cardiovascular Mortality", ylab="")
ts.plot(tempr, main="Temperature",  ylab="")
ts.plot(part, main="Particulates", ylab="")
#dev.new()  

#pairs(cbind( Mortality=cmort, Temperature=tempr, Particulates=part ))

#Regression
temp  = tempr-mean(tempr)  # center temperature    
temp2 = temp^2             # square it  
trend = time(cmort)        # time

fit = lm(cmort~ trend + temp + temp2 + part, na.action=NULL)
            
acf(resid(fit))


summary(fit)       # regression results
summary(aov(fit))  # ANOVA table   (compare to next line)
summary(aov(lm(cmort~cbind(trend, temp, temp2, part)))) # Table 2.1

num = length(cmort)                                     # sample size
AIC(fit)/num - log(2*pi)                                # AIC 
BIC(fit)/num - log(2*pi)                                # BIC 
# AIC(fit, k=log(num))/num - log(2*pi)                  # BIC (alt method)    
(AICc = log(sum(resid(fit)^2)/num) + (num+5)/(num-5-2)) # AICc
 


############################
#Section 2.2
#DIFFERENCING

#Differencing global temp
plot(globtemp, type="o")
acf(gtemp, 48, main="") 

par(mfrow=c(2,1))
ts.plot(diff(globtemp), type="o")
mean(diff(globtemp))     # drift estimate = .008
acf(diff(gtemp), 48, main="")


#Differencing chicken prices
fit = lm(chicken~time(chicken), na.action=NULL) # regress chicken on time
par(mfrow=c(2,1))
ts.plot(resid(fit), main="detrended")
ts.plot(diff(chicken), main="first difference")

 
par(mfrow=c(3,1))     # plot ACFs
acf(chicken, 48, main="chicken")
acf(resid(fit), 48, main="detrended")
acf(diff(chicken), 48, main="first difference")
#acf(1:100)


##############################
#VARIANCE-STABILIZATION

#Glacial varve thickness
par(mfrow=c(2,1))
ts.plot(varve, main="varve", ylab="")
ts.plot(log(varve), main="log(varve)", ylab="" )
 

#DOW JONES
dow=read.csv("C:\\Users\\cahoyd.UHD\\Documents\\Stat5307\\ShumwayStoffer\\F18\\Slides\\R Codes\\DowJones08302017b.csv", header = TRUE) 
dow2=dow[,2] #try 3
par(mfrow=c(2,2))
plot(dow2,type="l")
plot(log(dow2),type="l")
plot(diff(log(dow2)),type="l")
acf(diff(log(dow2)))

#OIL and GAS Data
par(mfrow=c(2,2))
plot(oil)
plot(log(oil))
plot(diff(log(oil)))
acf(diff(log(oil)))

par(mfrow=c(2,2))
plot(gas)
plot(log(gas))
plot(diff(log(gas)))
acf(diff(log(gas)))



######################################
#Discovering signal from noise
set.seed(1000)  # so you can reproduce these results
x = 2*cos(2*pi*1:500/50 + .6*pi) + rnorm(500,0,5)
z1 = cos(2*pi*1:500/50)  
z2 = sin(2*pi*1:500/50)
summary(fit <- lm(x~0+z1+z2))  # zero to exclude the intercept

par(mfrow=c(2,1))
ts.plot(x)
ts.plot(x, col=8, ylab=expression(hat(x)))
lines(fitted(fit), col=2) 
#dev.off()
 
acf(resid(fit))

############################
#Section 2.3: Smoothing

#Eg 2.11
#dev.new(height=4)
wgts = c(.5, rep(1,11), .5)/12
soif = filter(soi, sides=2, filter=wgts)
ts.plot(soi)
lines(soif, lwd=2, col=4)

#checking stationarity of residuals
acf2(soi-soif) 
 

#Eg 2.12
#dev.new(height=4)
ts.plot(soi)
trend=ksmooth(time(soi), soi, "normal", bandwidth=3.5)
lines(trend, lwd=2, col=4)
#par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert

#checking stationarity of residuals 
acf2(soi-trend$y)
 


#Eg 2.13
dev.off()
ts.plot(soi)
lines(lowess(soi, f=.05), lwd=2, col=4) # El Nino cycle
lines(lowess(soi,f=.5), lty=2, lwd=2, col=2) # trend (with default span f=2/3)



#Eg 2.14
ts.plot(soi)
lines(smooth.spline(time(soi), soi, spar=.5), lwd=2, col=4)
lines(smooth.spline(time(soi), soi, spar= 1), lty=2, lwd=2, col=2)

#Check stationarity of residuals for eg 2.13 & 2.14.



################
#SIGNAL DECOMPOSITION IN R
#install.packages("fpp")
library(fpp)
data(ausbeer)
head(ausbeer)
#Total quarterly beer production 
#in Australia (in megalitres) from 1956:Q1 to 2008:Q3
plot(ausbeer) 
#IMPORTANT FIRST STEP
#ausbeer=   ts(data,freq = ???) 
ausbeerdecomp=decompose(ausbeer, type="additive")
#ausbeerdecomp
#?decompose
plot(ausbeerdecomp)
ausbeerdecomp$figure # seasonal effects for each quarter

 
#More EDA
#install.packages("forecast")
library(forecast)
seasonplot(ausbeer)
monthplot(ausbeer)
boxplot(ausbeer ~ cycle(ausbeer))

#Another approach vis LOcal regrESSion
#Decompose a time series into seasonal, trend and irregular components using loess, acronym STL.
plot(stl(ausbeer, "periodic"))


#Johnson and Johnson Data
#J&J quarterly earnings
#load("tsa3.rda") # SEE THE FOOTNOTE
plot(jj, type="o", ylab="Quarterly Earnings per Share")
jjdecomp=decompose(jj, type="multiplicative")
jjdecomp$figure # seasonal effects for each quarter
plot(jjdecomp)

#EDA
seasonplot(jj)
monthplot(jj)
boxplot(jj ~ cycle(jj))

 

##################
#EXTRA
plot(decompose(rec))

 






