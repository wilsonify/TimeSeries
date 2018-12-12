# Exhibit 1.1 (CaC)
library(TSA)
#win.graph(width=4.875, height=2.5,pointsize=8)
data(larain)
plot(larain,ylab='Inches',xlab='Year',type='o', lwd=2, frame.plot=FALSE)
 
# Exhibit 1.3(CaC)
#win.graph(width=4.875, height=2.5,pointsize=8)
data(color)
plot(color,ylab='Color Property',xlab='Batch',type='o',lwd=2, frame.plot=FALSE)
 
# Exhibit 1.7(CaC)
#win.graph(width=4.875, height=2.5,pointsize=8)
data(tempdub)
plot(tempdub,ylab='Temperature',type='o')
 

 
#from SaS
library(astsa)
#J&J quarterly earnings
#load("tsa3.rda") # SEE THE FOOTNOTE
plot(jj, type="o", ylab="Quarterly Earnings per Share")

#Global temp
plot(gtemp, type="o", ylab="Global Temperature Deviations")


#El Nino and Fish pop'n
par(mfrow = c(2,1))  # set up the graphics
ts.plot(soi, ylab="", main="Southern Oscillation Index")
ts.plot(rec, ylab="", main="Recruitment") 

#White Noise
e=rnorm(100)
plot(e,type='l',ylab="White Noise",frame.plot=FALSE)
points(e)

# Exhibit 2.1 (CaC)
library(TSA)
#win.graph(width=4.875, height=2.5,pointsize=8)
# rwalk contains a simulated random walk without drift
data(rwalk)
plot(rwalk,type='o',ylab='Random Walk', frame.plot=FALSE)

# R code for simulating a random walk with, say 60, iid standard normal errors
n=60
set.seed(12345) # intialize the random number so that the simulation can be 
# reproducible.
rw=cumsum(rnorm(n))
sim.random.walk=ts(rw,freq=1,start=1)
plot(sim.random.walk,type='o',ylab='Another Random Walk')

#Random walk with or without drift
set.seed(154) # so you can reproduce the results
w = rnorm(200,0,1); x = cumsum(w) # two commands in one line
wd = w +.2; xd = cumsum(wd)
plot.ts(xd, ylim=c(-5,55), main="random walk with/out drift", frame.plot=FALSE)
lines(x); lines(.2*(1:200), lty="dashed")



#Illustrations of  sine and cosine curves
t=1:100; thet=2*pi*(t-1)/100
plot(cos(thet),typ="l", main="Cosine and sine curves", frame.plot=FALSE)
lines(sin(thet),typ="l", lty=2)


thet=2*pi*(t-1)/50 #halving the denominator doubles no. of cycles
plot(cos(thet),typ="l", main="Two cycles of sine and cosine waves", frame.plot=FALSE)
lines(sin(thet),typ="l", lty=2)




#Sinusoidal sine curve 
a=0.2 #amplitude
w=1/52 # frequency(one cycle/52 time points)
fi=0.7*pi #phase shift
t=0:150 #time t

par(mfrow=c(2,2))
mut=a*sin(2*pi*w*t +fi)
plot(t,mut,type="l",frame.plot=FALSE)
points(t,mut,type="p",pch=19)

yt=mut+rnorm(length(t) ) #N(0,1) noise
plot(t,yt,type="l",frame.plot=FALSE)
points(t,yt, pch=19)

yt2=mut+rnorm(length(t),0,4) #N(0,4^2)
plot(t,yt2,type="l",frame.plot=FALSE)
points(t,yt2, pch=19)

yt3=mut+rnorm(length(t),0,20) #N(20^2)
plot(t,yt3,type="l",frame.plot=FALSE)
points(t,yt3, pch=19)
dev.off()


#Random Cosine Curve (with random phase shift)
yt4=cos( 2*pi*( t/12 + runif(length(t)) ) )
plot(t,yt4,type="l",frame.plot=FALSE)
points(t,yt4, pch=19)


###################################
#Section 1.3-Section 1.5
library(astsa)
#AUTOCORRELATION FUNCTION:
#White Noise
n=150 # try increasing this
wn=rnorm(n)

par(mfrow=c(2,1))
plot(wn,type='l',ylab="White Noise",frame.plot=FALSE)
points(wn)
acf(wn, type="covariance",lag.max=10, plot=TRUE)

wn[1:5]

#Simple random walk
rw=cumsum(wn)
rw[1:4]
par(mfrow=c(2,1))
sim.random.walk=ts(rw,freq=1,start=1)
plot(sim.random.walk,type='o',ylab='Another Random Walk',frame.plot=FALSE)
acf(sim.random.walk, 50)$acf[1:10]


#Moving Average: MA
s.ma2=filter(wn,filter=rep(1/3, 3), method = "convolution",sides=1) #running average

par(mfrow=c(2,1))
sim.ma2=ts(s.ma2,freq=1,start=1)
plot(wn,type='o',ylab='Another Random Walk',frame.plot=FALSE)
lines(sim.ma2,type='o',ylab='Another Random Walk', col="red")
acf(na.omit(sim.ma2), 5)$acf
 
points(acf(na.omit(sim.ma2), 10)$acf[-1] )



#Eg 1.25
#From eg 1.5:
#par(mfrow = c(2,1))  # set up the graphics
ts.plot(soi, ylab="", main="Southern Oscillation Index")
#ts.plot(rec, ylab="", main="Recruitment") 

par(mfrow = c(2,1)) 
ts.plot(soi, ylab="", main="Southern Oscillation Index")
acf(soi,50)

#acf(soi, 6)
(r = round(acf(soi, 6, plot=FALSE)$acf, 3)) # first 7 sample acf values

#par(mfrow=c(1,2))
plot(lag(soi,-1), soi); legend('topleft', legend=r[1])
#plot(lag(soi,-6), soi); legend('topleft', legend=r[6])


#another acf example
data(sunspotz)
ts.plot(sunspotz)

eg2 = acf(sunspotz,36)
mat = cbind(0:36, eg2$acf)
mat 
colnames(mat) = c("lag", "acf")

#Another way with the lower and upper limit for alpha=0.05
plot(mat,type="h") 
abline(h=0)
abline(h=1.96/sqrt(459),col="blue")
abline(h=-1.96/sqrt(459),col="blue")

#p-value calculation
pv=2*pnorm( sqrt(459)*abs(mat[,2]),lower.tail = FALSE )
cbind(0:36, pv)




#Cross-correlation function
wn=rnorm(200)
xt=filter(wn,filter=rep(1, 2), method = "convolution",sides=1)
yt=filter(wn,filter=c(1, -1), method = "convolution",sides=1)
ccf(na.omit(xt),na.omit(yt), main="")

#compare
ccf(na.omit(yt),na.omit(xt), main="")

 
#Prediction using cross-correlation function
x=rnorm(100)
y=lag(x,-5) + rnorm(100)

ccf(x,y,ylab='CCorrF',type='correlation') #ccf(x,y,...)
#Note that negative lags indicate x leads y  

#ccf(y,x,ylab='CCovF',type='covariance') #ccf(y,x,...)
ccf(y,x,ylab='CCorrF',type='correlation') #ccf(y,x,...)
#Note that  positive  lags indicate x leads y 


#SOI and Recruitment Data
#  acf1 and ccf2 are astsa v1.7.7+ scripts
#  you can use acf and ccf instead
par(mfrow=c(3,1))
acf(soi, 48, main="Southern Oscillation Index")
acf(rec, 48, main="Recruitment")
soirec=ccf(rec, soi,48, main="SOI vs Recruitment") #ccf(x,y,...)

#Note: Peaks at negative lags indicate SOI(x)  leads Recruitment(y) 

#Lagged Regression: See Example 2.9: REC vs lagged-SOI for more details  
#lag-plots 
lag1.plot(soi, 12)
lag2.plot(soi,rec,8)

 



