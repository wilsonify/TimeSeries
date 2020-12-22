# Title     : Chapter 4
# Objective : ARMA Models
# Created by: thom
# Created on: 12/22/20


#Example 4.2

par(mfrow = c(2, 1))
tsplot(sarima.sim(ar = .9, n = 100), ylab = "x", col = 4, main = expression(AR(1) ~ ~~phi == +.9))
tsplot(sarima.sim(ar = -.9, n = 100), ylab = "x", col = 4, main = expression(AR(1) ~ ~~phi == -.9))

#Example 4.3

psi <- ARMAtoMA(ar = c(1.5, -.75), ma = 0, 50)
par(mfrow = c(2, 1), mar = c(2, 2.5, 1, 0) + .5, mgp = c(1.5, .6, 0), cex.main = 1.1)
plot(psi, xaxp = c(0, 144, 12), type = "n", col = 4, ylab = expression(psi - weights), main = expression(AR(2) ~ ~~phi[1] == 1.5 ~ ~~phi[2] == -.75))
abline(v = seq(0, 48, by = 12), h = seq(-.5, 1.5, .5), col = gray(.9))
lines(psi, type = "o", col = 4)
set.seed(8675309)
simulation <- arima.sim(list(order = c(2, 0, 0), ar = c(1.5, -.75)), n = 144)
plot(simulation, xaxp = c(0, 144, 12), type = "n", ylab = expression(X[~t]))
abline(v = seq(0, 144, by = 12), h = c(-5, 0, 5), col = gray(.9))
lines(simulation, col = 4)

#Examples 4.5

par(mfrow = c(2, 1))
tsplot(sarima.sim(ma = .9, n = 100), col = 4, ylab = "x", main = expression(MA(1) ~ ~~theta == +.9))
tsplot(sarima.sim(ma = -.9, n = 100), col = 4, ylab = "x", main = expression(MA(1) ~ ~~theta == -.9))

#Example 4.10

set.seed(8675309) # Jenny, I got your number
x <- rnorm(150, mean = 5) # generate iid N(5,1)s
arima(x, order = c(1, 0, 1))  # estimation

#Example 4.11

AR <- c(1, -.3, -.4) # original AR coefs on the left
polyroot(AR)
MA <- c(1, .5) # original MA coefs on the right
polyroot(MA)

#Example 4.12

round(ARMAtoMA(ar = .8, ma = -.5, 10), 2) # first 10 psi-weights
round(ARMAtoAR(ar = .8, ma = -.5, 10), 2) # first 10 pi-weights
ARMAtoMA(ar = 1, ma = 0, 20)

#Example 4.15

ACF <- ARMAacf(ar = c(1.5, -.75), ma = 0, 50)
plot(ACF, type = "h", xlab = "lag", panel.first = Grid())
abline(h = 0)

#Example 4.18

ACF <- ARMAacf(ar = c(1.5, -.75), ma = 0, 24)[-1]
PACF <- ARMAacf(ar = c(1.5, -.75), ma = 0, 24, pacf = TRUE)
par(mfrow = 1:2)
tsplot(ACF, type = "h", xlab = "lag", ylim = c(-.8, 1))
abline(h = 0)
tsplot(PACF, type = "h", xlab = "lag", ylim = c(-.8, 1))
abline(h = 0)

#Example 4.21

acf2(rec, 48) # will produce values and a graphic
(regr <- ar.ols(rec, order = 2, demean = FALSE, intercept = TRUE))
regr$asy.se.coef # standard errors of the estimates

#Example 4.24

rec.yw <- ar.yw(rec, order = 2)
rec.yw$x.mean # mean estimate
rec.yw$ar # phi parameter estimates
sqrt(diag(rec.yw$asy.var.coef)) # their standard errors
rec.yw$var.pred # error variance estimate

#Example 4.25

set.seed(1)
ma1 <- sarima.sim(ma = 0.9, n = 50)
acf1(ma1, plot = FALSE)[1]

#Example 4.27

tsplot(diff(log(varve)), col = 4, ylab = expression(nabla ~ log ~ X[~t]), main = "Transformed Glacial Varves")
acf2(diff(log(varve)))
#
x <- diff(log(varve)) # data
r <- acf1(x, 1, plot = FALSE) # acf(1)
c(0) -> w -> z -> Sc -> Sz -> Szw -> para # initialize
num <- length(x) # = 633
## Estimation
para[1] <- (1 - sqrt(1 - 4 * (r^2))) / (2 * r) # MME
niter <- 12
for (j in 1:niter) {
  for (i in 2:num) { w[i] <- x[i] - para[j] * w[i - 1]
    z[i] <- w[i - 1] - para[j] * z[i - 1]
  }
  Sc[j] <- sum(w^2)
  Sz[j] <- sum(z^2)
  Szw[j] <- sum(z * w)
  para[j + 1] <- para[j] + Szw[j] / Sz[j]
}
# Results
cbind(iteration = 1:niter - 1, thetahat = para[1:niter], Sc, Sz)
## Plot conditional SS
c(0) -> w -> cSS
th <- -seq(.3, .94, .01)
for (p in 1:length(th)) {
  for (i in 2:num) { w[i] <- x[i] - th[p] * w[i - 1]
  }
  cSS[p] <- sum(w^2)
}
tsplot(th, cSS, ylab = expression(S[c](theta)), xlab = expression(theta))
abline(v = para[1:12], lty = 2, col = 4) # add previous results to plot
points(para[1:12], Sc[1:12], pch = 16, col = 4)

#Example 4.28

sarima(diff(log(varve)), p = 0, d = 0, q = 1, no.constant = TRUE)

#Example 4.31

sarima(rec, p = 2, d = 0, q = 0) # fit the model
sarima.for(rec, n.ahead = 24, p = 2, d = 0, q = 0)
abline(h = 61.8585, col = 4) # display estimated mean
