# Title     : Chapter 3
# Objective : Time Series Regression and EDA
# Created by: thom
# Created on: 12/22/20


#Example 3.1

summary(fit <- lm(salmon ~ time(salmon), na.action = NULL))
tsplot(salmon, col = 4, ylab = "USD per KG", main = "Salmon Export Price")
abline(fit)

#Example 3.5

culer <- c(rgb(.66, .12, .85), rgb(.12, .66, .85), rgb(.85, .30, .12))
par(mfrow = c(3, 1))
tsplot(cmort, main = "Cardiovascular Mortality", col = culer[1], type = "o", pch = 19, ylab = "")
tsplot(tempr, main = "Temperature", col = culer[2], type = "o", pch = 19, ylab = "")
tsplot(part, main = "Particulates", col = culer[3], type = "o", pch = 19, ylab = "")
##
tsplot(cmort, main = "", ylab = "", ylim = c(20, 130), col = culer[1])
lines(tempr, col = culer[2])
lines(part, col = culer[3])
legend("topright", legend = c("Mortality", "Temperature", "Pollution"), lty = 1, lwd = 2, col = culer, bg = "white")
##
panel.cor <- function(x, y) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), 2)
  text(0.5, 0.5, r, cex = 1.75)
}

pairs(cbind(Mortality = cmort, Temperature = tempr, Particulates = part), col = "dodgerblue3", lower.panel = panel.cor)
##
par(mfrow = 2:1)
plot(tempr, tempr^2) # collinear
cor(tempr, tempr^2)
temp <- tempr - mean(tempr)
plot(temp, temp^2) # not collinear
cor(temp, temp^2)
##
temp <- tempr - mean(tempr) # center temperature
temp2 <- temp^2
trend <- time(cmort) # time
fit <- lm(cmort ~ trend + temp + temp2 + part, na.action = NULL)
summary(fit) # regression results
summary(aov(fit)) # ANOVA table (compare to next line)
summary(aov(lm(cmort ~ cbind(trend, temp, temp2, part)))) # Table 3.1
num <- length(cmort) # sample size
AIC(fit) / num - log(2 * pi) # AIC
BIC(fit) / num - log(2 * pi) # BIC

#Example 3.6

fish <- ts.intersect(rec, soiL6 = lag(soi, -6))
summary(fit1 <- lm(rec ~ soiL6, data = fish, na.action = NULL))
tsplot(resid(fit1), col = 4) # residuals
##
library(dynlm)
summary(fit2 <- dynlm(rec ~ L(soi, 6)))

#Example 3.10

fit <- lm(salmon ~ time(salmon), na.action = NULL) # the regression
par(mfrow = c(2, 1)) # plot transformed data
tsplot(resid(fit), main = "detrended salmon price")
tsplot(diff(salmon), main = "differenced salmon price")
par(mfrow = c(2, 1)) # plot their ACFs
acf1(resid(fit), 48, main = "detrended salmon price")
acf1(diff(salmon), 48, main = "differenced salmon price")

#Example 3.11

par(mfrow = c(2, 1))
tsplot(diff(gtemp_land), col = 4, main = "differenced global tmeperature")
mean(diff(gtemp_land)) # drift since 1880
# [1] 0.0143
acf1(diff(gtemp_land))
mean(window(diff(gtemp_land), start = 1980)) # drift since 1980
# [1] 0.0329

#Example 3.12

layout(matrix(1:4, 2), widths = c(2.5, 1))
par(mgp = c(1.6, .6, 0), mar = c(2, 2, .5, 0) + .5)
tsplot(varve, main = "", ylab = "", col = 4, margin = 0)
mtext("varve", side = 3, line = .5, cex = 1.2, font = 2, adj = 0)
tsplot(log(varve), main = "", ylab = "", col = 4, margin = 0)
mtext("log(varve)", side = 3, line = .5, cex = 1.2, font = 2, adj = 0)
qqnorm(varve, main = "", col = 4); qqline(varve, col = 2, lwd = 2)
qqnorm(log(varve), main = "", col = 4); qqline(log(varve), col = 2, lwd = 2)

#Example 3.13

lag1.plot(soi, 12, col = "dodgerblue3") # Figure 3.10
lag2.plot(soi, rec, 8, col = "dodgerblue3") # Figure 3.11

#Example 3.14

library(zoo)   # zoo allows easy use of the variable names
dummy <- ifelse(soi < 0, 0, 1)
fish <- as.zoo(ts.intersect(rec, soiL6 = lag(soi, -6), dL6 = lag(dummy, -6)))
summary(fit <- lm(rec ~ soiL6 * dL6, data = fish, na.action = NULL))
plot(fish$soiL6, fish$rec, panel.first = Grid(), col = 'dodgerblue3')
points(fish$soiL6, fitted(fit), pch = 3, col = 6)
lines(lowess(fish$soiL6, fish$rec), col = 4, lwd = 2)
tsplot(resid(fit))    # not shown
acf1(resid(fit))      # and obviously not noise

#Example 3.15

set.seed(90210) # so you can reproduce these results
x <- 2 * cos(2 * pi * 1:500 / 50 + .6 * pi) + rnorm(500, 0, 5)
z1 <- cos(2 * pi * 1:500 / 50)
z2 <- sin(2 * pi * 1:500 / 50)
summary(fit <- lm(x ~ 0 + z1 + z2)) # zero to exclude intercept
par(mfrow = c(2, 1))
tsplot(x, col = 4)
tsplot(x, ylab = expression(hat(x)), col = rgb(0, 0, 1, .5))
lines(fitted(fit), col = 2, lwd = 2)

#Example 3.16

w <- c(.5, rep(1, 11), .5) / 12
soif <- filter(soi, sides = 2, filter = w)
tsplot(soi, col = rgb(.5, .6, .85, .9), ylim = c(-1, 1.15))
lines(soif, lwd = 2, col = 4)
# insert
par(fig = c(.65, 1, .75, 1), new = TRUE)
w1 <- c(rep(0, 20), w, rep(0, 20))
plot(w1, type = "l", ylim = c(-.02, .1), xaxt = "n", yaxt = "n", ann = FALSE)

#Example 3.17

tsplot(soi, col = rgb(0.5, 0.6, 0.85, .9), ylim = c(-1, 1.15))
lines(ksmooth(time(soi), soi, "normal", bandwidth = 1), lwd = 2, col = 4)
# insert
par(fig = c(.65, 1, .75, 1), new = TRUE)
gauss <- function(x) { 1 / sqrt(2 * pi) * exp(-(x^2) / 2) }
curve(gauss(x), -3, 3, xaxt = "n", yaxt = "n", ann = FALSE, col = 4)

#
SOI <- ts(soi, freq = 1)
tsplot(SOI) # the time scale matters (not shown)
lines(ksmooth(time(SOI), SOI, "normal", bandwidth = 12), lwd = 2, col = 4)

#Example 3.18

tsplot(soi, col = rgb(0.5, 0.6, 0.85, .9))
lines(lowess(soi, f = .05), lwd = 2, col = 4) # El Niño cycle
# lines(lowess(soi), lty=2, lwd=2, col=2) # trend (with default span)
##-- trend with CIs using loess --##
lo <- predict(loess(soi ~ time(soi)), se = TRUE)
trnd <- ts(lo$fit, start = 1950, freq = 12) # put back ts attributes
lines(trnd, col = 6, lwd = 2)
L <- trnd - qt(0.975, lo$df) * lo$se
U <- trnd + qt(0.975, lo$df) * lo$se
xx <- c(time(soi), rev(time(soi)))
yy <- c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(.6, alpha = .4))

#Example 3.19

plot(tempr, cmort, xlab = "Temperature", ylab = "Mortality", col = 'dodgerblue3', panel.first = Grid())
lines(lowess(tempr, cmort), col = 4, lwd = 2)

#Example 3.20

x <- window(hor, start = 2002)
plot(decompose(x)) # not shown
plot(stl(x, s.window = "per")) # not shown
plot(stl(x, s.window = 15))  # nicer version below
##
culer <- c('cyan4', 4, 2, 6)
x <- window(hor, start = 2002)
par(mfrow = c(4, 1), cex.main = 1)
out <- stl(x, s.window = 15)$time.series
tsplot(x, main = 'Hawaiian Occupancy Rate', ylab = '% rooms', col = gray(.7))
text(x, labels = 1:4, col = culer, cex = 1.25)
tsplot(out[, 1], main = "Seasonal", ylab = '% rooms', col = gray(.7))
text(out[, 1], labels = 1:4, col = culer, cex = 1.25)
tsplot(out[, 2], main = "Trend", ylab = '% rooms', col = gray(.7))
text(out[, 2], labels = 1:4, col = culer, cex = 1.25)
tsplot(out[, 3], main = "Noise", ylab = '% rooms', col = gray(.7))
text(out[, 3], labels = 1:4, col = culer, cex = 1.25)
