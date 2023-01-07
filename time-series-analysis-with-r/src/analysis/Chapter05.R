# Title     : Chapter 5
# Objective : ARIMA Models
# Created by: thom
# Created on: 12/22/20


#Example 5.2

sarima(diff(log(varve)), p = 0, d = 0, q = 1, no.constant = TRUE)
# equivalently
sarima(log(varve), p = 0, d = 1, q = 1, no.constant = TRUE)

#Example 5.3

ARMAtoMA(ar = 1, ma = 0, 20) # psi-weights for rw

#Example 5.4

round(ARMAtoMA(ar = c(1.9, -.9), ma = 0, 60), 1) # ain't no AR(2)
#
set.seed(2001)
x <- sarima.sim(ar = .9, d = 1, n = 150)
y <- window(x, start = 1, end = 100)
sarima.for(y, n.ahead = 50, p = 1, d = 1, q = 0, plot.all = TRUE)
text(85, 255, "PAST"); text(115, 255, "FUTURE")
abline(v = 100, lty = 2, col = 4)
lines(x)

#Example 5.5

set.seed(666)
x <- sarima.sim(ma = -0.8, d = 1, n = 100)
(x.ima <- HoltWinters(x, beta = FALSE, gamma = FALSE))
plot(x.ima, main = "EWMA")

#Example 5.6

##-- Figure 5.3 --##
layout(1:2, heights = 2:1)
tsplot(gnp, col = 4)
acf1(gnp, main = "")
##-- Figure 5.4 --##
tsplot(diff(log(gnp)), ylab = "GNP Growth Rate", col = 4)
abline(mean(diff(log(gnp))), col = 6)
##-- Figure 5.5 --##
acf2(diff(log(gnp)), main = "")
##
sarima(diff(log(gnp)), 0, 0, 2) # MA(2) on growth rate
sarima(diff(log(gnp)), 1, 0, 0) # AR(1) on growth rate
#
round(ARMAtoMA(ar = .35, ma = 0, 10), 3) # print psi-weights

#Example 5.7

sarima(diff(log(gnp)), 0, 0, 2) # MA(2) fit with diagnostics

#Example 5.8

sarima(log(varve), 0, 1, 1, no.constant = TRUE) # ARIMA(0,1,1)
sarima(log(varve), 1, 1, 1, no.constant = TRUE) # ARIMA(1,1,1)

#Example 5.9

uspop <- c(75.995, 91.972, 105.711, 123.203, 131.669, 150.697, 179.323, 203.212, 226.505, 249.633, 281.422, 308.745)
uspop <- ts(uspop, start = 1900, freq = .1)
t <- time(uspop) - 1955
reg <- lm(uspop ~ t +
  I(t^2) +
  I(t^3) +
  I(t^4) +
  I(t^5) +
  I(t^6) +
  I(t^7) +
  I(t^8))
b <- as.vector(reg$coef)

g <- function(t) { b[1] +
  b[2] * (t - 1955) +
  b[3] * (t - 1955)^2 +
  b[4] * (t - 1955)^3 +
  b[5] * (t - 1955)^4 +
  b[6] * (t - 1955)^5 +
  b[7] * (t - 1955)^6 +
  b[8] * (t - 1955)^7 +
  b[9] * (t - 1955)^8
}

par(mar = c(2, 2.5, .5, 0) + .5, mgp = c(1.6, .6, 0))
curve(g, 1900, 2024, ylab = "Population", xlab = "Year", main = "U.S. Population by Official Census", panel.first = Grid(), cex.main = 1, font.main = 1, col = 4)
abline(v = seq(1910, 2020, by = 20), lty = 1, col = gray(.9))
points(time(uspop), uspop, pch = 21, bg = rainbow(12), cex = 1.25)
mtext(expression("" %*% 10^6), side = 2, line = 1.5, adj = .95)
axis(1, seq(1910, 2020, by = 20), labels = TRUE)

#Example 5.10

sarima(diff(log(gnp)), 1, 0, 0) # AR(1)
sarima(diff(log(gnp)), 0, 0, 2) # MA(2)

#Example 5.11

set.seed(666)
phi <- c(rep(0, 11), .9)
sAR <- ts(arima.sim(list(order = c(12, 0, 0), ar = phi), n = 37), freq = 12) + 50
layout(matrix(c(1, 2, 1, 3), nc = 2), heights = c(1.5, 1))
par(mar = c(2.5, 2.5, 2, 1), mgp = c(1.6, .6, 0))
plot(sAR, xaxt = "n", col = gray(.6), main = "seasonal AR(1)", xlab = "YEAR", type = "c", ylim = c(45, 54))
abline(v = 1:4, lty = 2, col = gray(.6))
axis(1, 1:4); box()
abline(h = seq(46, 54, by = 2), col = gray(.9))
Months <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
points(sAR, pch = Months, cex = 1.35, font = 4, col = 1:4)
ACF <- ARMAacf(ar = phi, ma = 0, 100)[-1]
PACF <- ARMAacf(ar = phi, ma = 0, 100, pacf = TRUE)
LAG <- 1:100 / 12
plot(LAG, ACF, type = "h", xlab = "LAG", ylim = c(-.1, 1), axes = FALSE)
segments(0, 0, 0, 1)
axis(1, seq(0, 8, by = 1)); axis(2); box(); abline(h = 0)
plot(LAG, PACF, type = "h", xlab = "LAG", ylim = c(-.1, 1), axes = FALSE)
axis(1, seq(0, 8, by = 1)); axis(2); box(); abline(h = 0)

#Example 5.12

##-- Figure 5.10 --##
phi <- c(rep(0, 11), .8)
ACF <- ARMAacf(ar = phi, ma = -.5, 50)[-1]
PACF <- ARMAacf(ar = phi, ma = -.5, 50, pacf = TRUE)
LAG <- 1:50 / 12
par(mfrow = c(1, 2))
plot(LAG, ACF, type = "h", ylim = c(-.4, .8), panel.first = Grid())
abline(h = 0)
plot(LAG, PACF, type = "h", ylim = c(-.4, .8), panel.first = Grid())
abline(h = 0)
##-- birth series --##
tsplot(birth) # monthly number of births in US
acf2(diff(birth)) # P/ACF of the differenced birth rate

##-- seasonal persistence --##
x <- window(hor, start = 2002)
par(mfrow = c(2, 1))
tsplot(x, main = "Hawaiian Quarterly Occupancy Rate", ylab = " % rooms", ylim = c(62, 86), col = gray(.7))
text(x, labels = 1:4, col = c(3, 4, 2, 6), cex = .8)
Qx <- stl(x, 15)$time.series[, 1]
tsplot(Qx, main = "Seasonal Component", ylab = " % rooms", ylim = c(-4.7, 4.7), col = gray(.7))
text(Qx, labels = 1:4, col = c(3, 4, 2, 6), cex = .8)

#Example 5.15

par(mfrow = c(2, 1))
tsplot(cardox, col = 4, ylab = expression(CO[2]))
title("Monthly Carbon Dioxide Readings - Mauna Loa Observatory ", cex.main = 1)
tsplot(diff(diff(cardox, 12)), col = 4,
       ylab = expression(nabla ~ nabla[12] ~ CO[2]))

acf2(diff(diff(cardox, 12)))

sarima(cardox, p = 0, d = 1, q = 1, P = 0, D = 1, Q = 1, S = 12)
sarima(cardox, p = 1, d = 1, q = 1, P = 0, D = 1, Q = 1, S = 12)

sarima.for(cardox, 60, 1, 1, 1, 0, 1, 1, 12)
abline(v = 2018.9, lty = 6)
##-- for comparison --##
sarima.for(cardox, 60, 0, 1, 1, 0, 1, 1, 12) # not shown

#Example 5.16

trend <- time(cmort); temp <- tempr - mean(tempr); temp2 <- temp^2
fit <- lm(cmort ~ trend + temp + temp2 + part, na.action = NULL)
acf2(resid(fit), 52) # implies AR2
sarima(cmort, 2, 0, 0, xreg = cbind(trend, temp, temp2, part))

#Example 5.17

library(zoo)
lag2.plot(Hare, Lynx, 5) # lead-lag relationship
pp <- as.zoo(ts.intersect(Lynx, HareL1 = lag(Hare, -1)))
summary(reg <- lm(pp$Lynx ~ pp$HareL1)) # results not displayed
acf2(resid(reg)) # in Figure 5.11
(reg2 <- sarima(pp$Lynx, 2, 0, 0, xreg = pp$HareL1))
prd <- Lynx - resid(reg2$fit) # prediction (resid = obs - pred)
prde <- sqrt(reg2$fit$sigma2) # prediction error
tsplot(prd, lwd = 2, col = rgb(0, 0, .9, .5), ylim = c(-20, 90), ylab = "Lynx")
points(Lynx, pch = 16, col = rgb(.8, .3, 0))
x <- time(Lynx)[-1]
xx <- c(x, rev(x))
yy <- c(prd - 2 * prde, rev(prd + 2 * prde))
polygon(xx, yy, border = 8, col = rgb(.4, .5, .6, .15))
mtext(expression("" %*% 10^3), side = 2, line = 1.5, adj = .975)
legend("topright", legend = c("Predicted", "Observed"), lty = c(1, NA), lwd = 2, pch = c(NA, 16), col = c(4, rgb(.8, .3, 0)), cex = .9)
