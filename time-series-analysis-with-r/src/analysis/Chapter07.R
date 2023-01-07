# Title     : Chapter 7
# Objective : Spectral Estimation
# Created by: thom
# Created on: 12/22/20


DFTs

(dft <- fft(1:4) / sqrt(4))
(idft <- fft(dft, inverse = TRUE) / sqrt(4))

#Example 7.4

par(mfrow = c(2, 1)) # raw periodogram
mvspec(soi, col = rgb(.05, .6, .75), lwd = 2)
rect(1 / 7, -1e5, 1 / 3, 1e5, density = NA, col = gray(.5, .2))
abline(v = 1 / 4, lty = 2, col = "dodgerblue")
mtext("1/4", side = 1, line = 0, at = .25, cex = .75)
mvspec(rec, col = rgb(.05, .6, .75), lwd = 2)
rect(1 / 7, -1e5, 1 / 3, 1e5, density = NA, col = gray(.5, .2))
abline(v = 1 / 4, lty = 2, col = "dodgerblue")
mtext("1/4", side = 1, line = 0, at = .25, cex = .75)

#  log redux
par(mfrow = c(2, 1)) # raw periodogram
mvspec(soi, col = rgb(.05, .6, .75), lwd = 2, log = 'y')
rect(1 / 7, 1e-5, 1 / 3, 1e5, density = NA, col = gray(.5, .2))
abline(v = 1 / 4, lty = 2, col = "dodgerblue")
mtext("1/4", side = 1, line = 0, at = .25, cex = .75)
mvspec(rec, col = rgb(.05, .6, .75), lwd = 2, log = 'y')
rect(1 / 7, 1e-5, 1 / 3, 1e5, density = NA, col = gray(.5, .2))
abline(v = 1 / 4, lty = 2, col = "dodgerblue")
mtext("1/4", side = 1, line = 0, at = .25, cex = .75)

Periodogram...Bad!

u = fft(rnorm(2^10)) # DFT of the data
z = Mod(u/2^5)^2 # periodogram
w = 0:511/1024 # frequencies
tsplot(w, z[1:512], col=rgb(.05, .6, .75), ylab="Periodogram", xlab="Frequency")
segments(0, 1, .5, 1, col=rgb(1, .25, 0), lwd=5) # actual spectrum
fz = filter(z, filter=rep(.01, 100), circular=TRUE)# smooth
lines(w, fz[1:512], col=rgb(0, .25, 1, .7), lwd=3) # plot the smooth
#
# Get the same thing with less code and n = 1000
u = mvspec(rnorm(1000), col=rgb(.05, .6, .75)) # periodogram
abline(h=1, col=2, lwd=5) # true spectrum
lines(u$freq, filter(u$spec, filter=rep(1, 101)/101, circular=TRUE), col=4, lwd=2) # add the smooth

#Example 7.5

par(mfrow=c(2, 1))
soi.ave = mvspec(soi, spans=9, col=rgb(.05, .6, .75), lwd=2)
rect(1/7, -1e5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)
rec.ave = mvspec(rec, spans=9, col=rgb(.05, .6, .75), lwd=2)
rect(1/7, -1e5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)

##-- redo on log scale with CIs --##
par(mfrow=c(2, 1))
soi.ave = mvspec(soi, spans=9, col=rgb(.05, .6, .75), lwd=2, log='yes')
rect(1/7, 1e-5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)
rec.ave = mvspec(rec, spans=9, col=rgb(.05, .6, .75), lwd=2, log='yes')
rect(1/7, 1e-5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)

#Example 7.6

y = ts(rev(1:100 %% 20), freq=20) # sawtooth signal
par(mfrow=2:1)
tsplot(1:100, y, ylab="sawtooth signal", col=4)
mvspec(y, main="", ylab="periodogram", col=rgb(.05, .6, .75), xlim=c(0, 7))

#Example 7.7

(dm = kernel("modified.daniell", c(3, 3))) # for a list
# the figure with both kernels
par(mfrow=1:2, mar=c(3, 3, 2, 1), mgp=c(1.6, .6, 0))
plot(kernel("modified.daniell", c(3, 3)), ylab=expression(h[~k]), cex.main=1, col=4, panel.first=Grid())
plot(kernel("modified.daniell", c(3, 3, 3)), ylab=expression(h[~k]), cex.main=1, col=4, panel.first=Grid())
#
par(mfrow=c(2, 1))
sois = mvspec(soi, spans=c(7, 7), taper=.1, col=rgb(.05, .6, .75), lwd=2)
rect(1/7, -1e5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)
recs = mvspec(rec, spans=c(7, 7), taper=.1, col=rgb(.05, .6, .75), lwd=2)
rect(1/7, -1e5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)
sois$Lh
sois$bandwidth
# to find the peaks
sois$details[1:45,]

##-- for the logs - not shown in the text --##
par(mfrow=c(2, 1))
sois = mvspec(soi, spans=c(7, 7), taper=.1, col=rgb(.05, .6, .75), lwd=2, log='yes')
rect(1/7, 1e-5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)
recs = mvspec(rec, spans=c(7, 7), taper=.1, col=rgb(.05, .6, .75), lwd=2, log='yes')
rect(1/7, 1e-5, 1/3, 1e5, density=NA, col=gray(.5, .2))
abline(v=.25, lty=2, col="dodgerblue")
mtext("1/4", side=1, line=0, at=.25, cex=.75)

Tapering

w = seq(-.04, .04, .0001); n=480; u=0
for (i in -4:4){
k = i/n
u = u + sin(n*pi*(w+k))^2 / sin(pi*(w+k))^2
}
fk = u/(9*480)
u=0; wp = w+1/n; wm = w-1/n
for (i in -4:4){
k  = i/n; wk = w+k; wpk = wp+k; wmk = wm+k
z  = complex(real=0, imag=2*pi*wk)
zp = complex(real=0, imag=2*pi*wpk)
zm = complex(real=0, imag=2*pi*wmk)
d  = exp(z)*(1-exp(z*n))/(1-exp(z))
dp = exp(zp)*(1-exp(zp*n))/(1-exp(zp))
dm = exp(zm)*(1-exp(zm*n))/(1-exp(zm))
D  = .5*d - .25*dm*exp(pi*w/n)-.25*dp*exp(-pi*w/n)
D2 = abs(D)^2
u  = u + D2
}
sfk =u/(480*9)
par(mfrow=c(1, 2))
plot(w, fk, type="l", ylab="", xlab="frequency", main="Without Tapering", yaxt="n")
mtext(expression("|"), side=1, line=-.20, at=c(-0.009375, .009375), cex=1.5, col=2)
segments(-4.5/480, -2, 4.5/480, -2, lty=1, lwd=3, col=2)
plot(w, sfk, type="l", ylab="", xlab="frequency", main="With Tapering", yaxt="n")
mtext(expression("|"), side=1, line=-.20, at=c(-0.009375, .009375), cex=1.5, col=2)
segments(-4.5/480, -.78, 4.5/480, -.78, lty=1, lwd=3, col=2)

#Example 7.8

par(mar=c(2.5, 2.5, 1, 1), mgp=c(1.5, .6, 0))
s0 = mvspec(soi, spans=c(7, 7), plot=FALSE)# no taper
s20 = mvspec(soi, spans=c(7, 7), taper=.2, plot=FALSE) # 20% taper
s50 = mvspec(soi, spans=c(7, 7), taper=.5, plot=FALSE) # full taper
plot(s0$freq[1:70], s0$spec[1:70], log="y", type="l", ylab="log-spectrum", xlab="frequency", panel.first=Grid())
lines(s20$freq[1:70], s10$spec[1:70], col=2)
lines(s50$freq[1:70], s50$spec[1:70], col=4)
text(.72, 0.04, "leakage", cex=.8)
arrows(.72, .035, .70, .011, length=0.05, angle=30)
abline(v=.25, lty=2, col=8)
mtext("1/4", side=1, line=0, at=.25, cex=.9)
legend("bottomleft", legend=c("no taper", "20% taper", "50% taper"), lty=1, col=c(1, 2, 4), bty="n")
par(fig = c(.7, 1, .7, 1), new = TRUE)
taper <- function(x) { .5*(1+cos(2*pi*x)) }
x <- seq(from = -.5, to = .5, by = 0.001)
plot(x, taper(x), type="l", lty=1, yaxt="n", xaxt="n", ann=FALSE)

#Example 7.10

spaic = spec.ar(soi, log="no", col="cyan4") # min AIC spec
abline(v=frequency(soi)*1/48, lty="dotted") # El Niño Cycle
(soi.ar = ar(soi, order.max=30)) # estimates and AICs
plot(1:30, soi.ar$aic[-1], type="o") # plot AICs
##
n = length(soi)
c()-> AIC -> BIC
for (k in 1:30){
sigma2 = ar(soi, order=k, aic=FALSE)$var.pred
BIC[k] = log(sigma2) + k*log(n)/n
AIC[k] = log(sigma2)+ (n+2*k)/n
IC = cbind(AIC, BIC+1)
ts.plot(IC, type="o", xlab="p", ylab="AIC / BIC")
Grid()
}

#Example 7.12

sr = mvspec(cbind(soi, rec), kernel("daniell", 9), plot=FALSE)
sr$df
(f = qf(.999, 2, sr$df-2))
(C = f/(18+f))
plot(sr, plot.type = "coh", ci.lty = 2, main="SOI & Recruitment")
abline(h = C)
