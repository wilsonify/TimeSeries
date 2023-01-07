# Title     : Chapter 1
# Objective : Time Series Elements
# Created by: thom
# Created on: 12/22/20


#Example 1.1

par(mfrow = 2:1)
tsplot(jj, ylab = "QEPS", type = "o", col = 4, main = "Johnson & Johnson Quarterly Earnings")
tsplot(log(jj), ylab = "log(QEPS)", type = "o", col = 4)

#Example 1.2

culer <- c(rgb(217, 77, 30, 160, max = 255), rgb(30, 170, 217, 160, max = 255))
tsplot(gtemp_land, col = culer[1], lwd = 2, type = "o", pch = 20,
       ylab = "Temperature Deviations", main = "Global Warming")
lines(gtemp_ocean, col = culer[2], lwd = 2, type = "o", pch = 20)
legend("topleft", col = culer, lty = 1, lwd = 2, pch = 20, legend = c("Land Surface", "Sea Surface"), bg = "white")

#Example 1.3

library(xts)   # install.packages("xts") if you don't have it already
djia_return <- diff(log(djia$Close))[-1]
par(mfrow = 2:1)
plot(djia$Close, col = 4)
plot(djia_return, col = 4)

tsplot(diff(log(gdp)), type = "o", col = 4) # using diff log
points(diff(gdp) / lag(gdp, -1), pch = "+", col = 2) # actual return

#Example 1.4

par(mfrow = c(2, 1))
tsplot(soi, ylab = "", xlab = "", main = "Southern Oscillation Index", col = 4)
text(1970, .91, "COOL", col = "cyan4")
text(1970, -.91, "WARM", col = "darkmagenta")
tsplot(rec, ylab = "", main = "Recruitment", col = 4)

#Example 1.5

culer <- c(rgb(.85, .30, .12, .6), rgb(.12, .67, .86, .6))
tsplot(Hare, col = culer[1], lwd = 2, type = "o", pch = 0,
       ylab = expression(Number ~ ~~("" %*% 1000)))
lines(Lynx, col = culer[2], lwd = 2, type = "o", pch = 2)
legend("topright", col = culer, lty = 1, lwd = 2, pch = c(0, 2), legend = c("Hare", "Lynx"), bty = "n")

#Example 1.6

par(mfrow = c(3, 1))
culer <- c(rgb(.12, .67, .85, .7), rgb(.67, .12, .85, .7))
u <- rep(c(rep(.6, 16), rep(-.6, 16)), 4) # stimulus signal
tsplot(fmri1[, 4], ylab = "BOLD", xlab = "", main = "Cortex", col = culer[1], ylim = c(-.6, .6), lwd = 2)
lines(fmri1[, 5], col = culer[2], lwd = 2)
lines(u, type = "s")
tsplot(fmri1[, 6], ylab = "BOLD", xlab = "", main = "Thalamus", col = culer[1], ylim = c(-.6, .6), lwd = 2)
lines(fmri1[, 7], col = culer[2], lwd = 2)
lines(u, type = "s")
tsplot(fmri1[, 8], ylab = "BOLD", xlab = "", main = "Cerebellum", col = culer[1], ylim = c(-.6, .6), lwd = 2)
lines(fmri1[, 9], col = culer[2], lwd = 2)
lines(u, type = "s")
mtext("Time (1 pt = 2 sec)", side = 1, line = 1.75)

#Example 1.7 - 1.8

par(mfrow = 2:1)
w <- rnorm(500) # 500 N(0,1) variates
v <- filter(w, sides = 2, filter = rep(1 / 3, 3)) # moving average
tsplot(w, col = 4, main = "white noise")
tsplot(v, ylim = c(-3, 3), col = 4, main = "white noise")

#Example 1.9

set.seed(90210)
w <- rnorm(250 + 50) # 50 extra to avoid startup problems
x <- filter(w, filter = c(1.5, -.75), method = "recursive")[-(1:50)]
tsplot(x, main = "autoregression", col = 4)

#Example 1.10

set.seed(314159265) # so you can reproduce the results
w <- rnorm(200)
x <- cumsum(w)
wd <- w + .3
xd <- cumsum(wd)
tsplot(xd, ylim = c(-2, 80), main = "random walk", ylab = "", col = 4)
abline(a = 0, b = .3, lty = 2, col = 4) # drift
lines(x, col = "darkred")
abline(h = 0, col = "darkred", lty = 2)

#Example 1.11

cs <- 2 * cos(2 * pi * (1:500) / 50 + .6 * pi)
w <- rnorm(500, 0, 1)
par(mfrow = c(3, 1), mar = c(3, 2, 2, 1), cex.main = 1.5)   # help(par) for info
tsplot(cs, ylab = "", main = expression(x[t] == 2 * cos(2 * pi * t / 50 + .6 * pi)))
tsplot(cs + w, ylab = "", main = expression(x[t] == 2 * cos(2 * pi * t / 50 + .6 * pi) + N(0, 1)))
tsplot(cs + 5 * w, ylab = "", main = expression(x[t] == 2 * cos(2 * pi * t / 50 + .6 * pi) + N(0, 25)))
