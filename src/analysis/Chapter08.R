# Title     : Chapter 8
# Objective : GARCH Models, Unit Root Testing, Long Memory and Fractional Differencing, State Space Models, Cross-Correlation Analysis and Prewhitening, Bootstrapping Autoregressive Models, Threshold Autoregressive Models
# Created by: thom
# Created on: 12/22/20


#Example 8.1

res <- resid(sarima(diff(log(gnp)), 1, 0, 0, details = FALSE)$fit)
acf2(res^2, 20)
#
library(fGarch)
gnpr <- diff(log(gnp))
summary(garchFit(~arma(1, 0) + garch(1, 0), data = gnpr))

#Example 8.2

library(xts)
djiar <- diff(log(djia$Close))[-1]
acf2(djiar) # exhibits some autocorrelation
u <- resid(sarima(djiar, 1, 0, 0, details = FALSE)$fit)
acf2(u^2)   # oozes autocorrelation
library(fGarch)
summary(djia.g <- garchFit(~arma(1, 0) + garch(1, 1), data = djiar, cond.dist = "std"))
plot(djia.g, which = 3)

#Example 8.3

lapply(c("xts", "fGarch"), library, char = TRUE) # load 2 packages in one line - amazing!
djiar <- diff(log(djia$Close))[-1]
summary(djia.ap <- garchFit(~arma(1, 0) + aparch(1, 1), data = djiar, cond.dist = "std"))
plot(djia.ap)   # to see all plot options

#Example 8.4

layout(1:2)
acf1(cumsum(rnorm(634)), 100, main = "Series: random walk")
acf1(log(varve), 100, ylim = c(-.1, 1))
#
library(tseries)
adf.test(log(varve), k = 0) # DF test
adf.test(log(varve))      # ADF test
pp.test(log(varve))       # PP test

#Example 8.5

d <- 0.3727893; p <- c(1)
for (k in 1:30) { p[k + 1] <- (k - d) * p[k] / (k + 1) }
tsplot(1:30, p[-1], ylab = expression(pi(d)), lwd = 2, xlab = "Index", type = "h", col = "dodgerblue3")
library(arfima)
summary(varve.fd <- arfima(log(varve), order = c(0, 0, 0)))
# residuals
innov <- resid(varve.fd)
tsplot(innov[[1]]) # not shown
par(mfrow = 2:1)
acf1(resid(sarima(log(varve), 1, 1, 1, details = FALSE)$fit), main = "ARIMA(1,1,1)")
acf1(innov[[1]], main = "Frac Diff")

#Example 8.8

u <- ssm(gtemp_land, A = 1, alpha = .01, phi = 1, sigw = .01, sigv = .1)
tsplot(gtemp_land, col = "dodgerblue3", type = "o", pch = 20, ylab = "Temperature Deviations")
lines(u$Xs, col = 6, lwd = 2)
xx <- c(time(u$Xs), rev(time(u$Xs)))
yy <- c(u$Xs - 2 * sqrt(u$Ps), rev(u$Xs + 2 * sqrt(u$Ps)))
polygon(xx, yy, border = 8, col = gray(.6, alpha = .25))

#Example 8.9

ccf2(cmort, part) # Figure 8.7
acf2(diff(cmort)) # Figure 8.8 implies AR(1)
u <- sarima(cmort, 1, 1, 0, no.constant = TRUE) # fits well
cmortw <- resid(u$fit)
phi <- as.vector(u$fit$coef) # -.5064
# filter particluates the same way
partf <- filter(diff(part), filter = c(1, -phi), sides = 1)
## -- now line up the series - this step is important --##
both <- ts.intersect(cmortw, partf) # line them up
Mw <- both[, 1] # cmort whitened
Pf <- both[, 2] # part filtered
ccf2(Mw, Pf) # Figure 8.9

--Section8.6--

# data
set.seed(101010)
e = rexp(150, rate=.5); u = runif(150, -1, 1); de = e*sign(u)
dex = 50 + arima.sim(n=100, list(ar=.95), innov=de, n.start=50)
layout(matrix(1:2, nrow=1), widths=c(5, 2))
tsplot(dex, col=4, ylab=expression(X[~t]))
# density - standard Laplace vs normal
f = function(x) { .5*dexp(abs(x), rate = 1/sqrt(2))}
curve(f, -5, 5, panel.first=Grid(), col=4, ylab="f(w)", xlab="w")
par(new=TRUE)
curve(dnorm, -5, 5, ylab="", xlab="", yaxt="no", xaxt="no", col=2)
#
fit = ar.yw(dex, order=1)
round(cbind(fit$x.mean, fit$ar, fit$var.pred), 2)
#
set.seed(111)
phi.yw = c()
for (i in 1:1000){
e  = rexp(150, rate=.5)
u  = runif(150, -1, 1)
de = e*sign(u)
x  = 50 + arima.sim(n=100, list(ar=.95), innov=de, n.start=50)
phi.yw[i] = ar.yw(x, order=1)$ar
}
#
set.seed(666) # not that 666
fit = ar.yw(dex, order=1) # assumes the data were retained
m = fit$x.mean # estimate of mean
phi = fit$ar # estimate of phi
nboot = 500 # number of bootstrap replicates
resids = fit$resid[-1] # the 99 residuals
x.star = dex # initialize x*
phi.star.yw = c()
# Bootstrap
for (i in 1:nboot) {
resid.star = sample(resids, replace=TRUE)
for (t in 1:99)
{
x.star[t+1] = m + phi*(x.star[t]-m)+ resid.star[t]
}
phi.star.yw[i]= ar.yw(x.star, order=1)$ar
}
# Picture
culer = rgb(0, .5, .5, .5)
hist(phi.star.yw, 15, main="", prob=TRUE, xlim=c(.65, 1.05),
ylim=c(0, 14), col=culer, xlab=expression(hat(phi)))
lines(density(phi.yw, bw=.02), lwd=2) # from previous simulation
u = seq(.75, 1.1, by=.001)# normal approximation
lines(u, dnorm(u, mean=.96, sd=.03), lty=2, lwd=2)
legend(.65, 14, legend=c("true distribution", "bootstrap distribution", "normal approximation"), bty="n", lty=c(1, 0, 2), lwd=c(2, 0, 2), col=1, pch=c(NA, 22, NA), pt.bg=c(NA, culer, NA), pt.cex=2.5)
# CIs
alf = .025 # 95% CI
quantile(phi.star.yw, probs = c(alf, 1-alf))
quantile(phi.yw, probs = c(alf, 1-alf))
n=100; phi = fit$ar
se = sqrt((1-phi)/n)
c(phi - qnorm(1-alf)*se, phi + qnorm(1-alf)*se)

#Example 8.10

tsplot(flu, type="c", ylab="Influenza Deaths per 10,000")
Months = c("J", "F", "M", "A", "M","J", "J", "A", "S", "O","N", "D")
culers = c(rgb(0, .4, .8), rgb(.8, 0, .4), rgb(0, .8, .4), rgb(.8, .4, 0))
points(flu, pch=Months, cex=.8, font=4, col=culers)
# Start analysis
dflu = diff(flu)
lag1.plot(dflu, corr=FALSE) # scatterplot with lowess fit
thrsh = .05 # threshold
Z = ts.intersect(dflu, lag(dflu, -1), lag(dflu, -2), lag(dflu, -3),
lag(dflu, -4))
ind1 = ifelse(Z[, 2] < thrsh, 1, NA) # indicator < thrsh
ind2 = ifelse(Z[, 2] < thrsh, NA, 1)# indicator >= thrsh
X1 = Z[, 1]*ind1
X2= Z[, 1]*ind2
summary(fit1 <- lm(X1~ Z[, 2:5])) # case 1
summary(fit2 <- lm(X2~ Z[, 2:5])) # case 2
D = cbind(rep(1, nrow(Z)), Z[, 2:5]) # design matrix
p1 = D %*% coef(fit1) # get predictions
p2 = D %*% coef(fit2)
prd = ifelse(Z[, 2] < thrsh, p1, p2)
# Figure 8.11
tsplot(prd, ylim=c(-.5, .5), ylab=expression(nabla~flu[~t]), lwd=2,
col=rgb(0, 0, .9, .5))
prde1 = sqrt(sum(resid(fit1)^2)/df.residual(fit1))
prde2 = sqrt(sum(resid(fit2)^2)/df.residual(fit2))
prde = ifelse(Z[, 2] < thrsh, prde1, prde2)
x = time(dflu)[-(1:4)]
x = c(x, rev(x))
yy = c(prd - 2*prde, rev(prd + 2*prde))
polygon(xx, yy, border=8, col=rgb(.4, .5, .6, .15))
abline(h=.05, col=4, lty=6)
points(dflu, pch=16, col="darkred")
#
par(mar=c(2.5, 2.5, 0, 0)+.5, mgp=c(1.6, .6, 0))
U = matrix(Z, ncol=5) # Z was created in the analysis above
culer = c(rgb(0, 1, 0, .4), rgb(0, 0, 1, .4))
culers = ifelse(U[, 2]<.05, culer[1], culer[2])
plot(U[, 2], U[, 1], panel.first=Grid(), pch=21, cex=1.1, bg=culers, xlab=expression(nabla~flu[~t-1]), ylab=expression(nabla~flu[~t]))
lines(lowess(U[, 2], U[, 1], f=2/3), col=6)
abline(v=.05, lty=2, col=4)
##
library(tsDyn) # load package - install it if you don"t have it
# vignette("tsDyn") # for package details
(u = setar(dflu, m=4, thDelay=0, th=.05)) # fit model and view results
(u = setar(dflu, m=4, thDelay=0))# let program fit threshold (=.036)
BIC(u); AIC(u)# if you want to try other models; m=3 works well too
plot(u) # graphics - ?plot.setar for information

