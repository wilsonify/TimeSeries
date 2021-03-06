# Title     : Chapter 2
# Objective : Time Series Regression and Exploratory Data Analysis
# Created by: thom
# Created on: 12/22/20


#Example 2.1

summary(fit <- lm(chicken ~ time(chicken))) # regress price on time
tsplot(chicken, ylab = "cents per pound", col = 4, lwd = 2)
abline(fit)           # add the fitted regression line to the plot

#Example 2.2

par(mfrow = c(3, 1))
tsplot(cmort, main = "Cardiovascular Mortality", ylab = "")
tsplot(tempr, main = "Temperature", ylab = "")
tsplot(part, main = "Particulates", ylab = "")

dev.new()
pairs(cbind(Mortality = cmort, Temperature = tempr, Particulates = part))

#  Regression
temp <- tempr - mean(tempr)  # center temperature
temp2 <- temp^2             # square it
trend <- time(cmort)        # time

fit <- lm(cmort ~ trend + temp + temp2 + part, na.action = NULL)

summary(fit)       # regression results
summary(aov(fit))  # ANOVA table   (compare to next line)
summary(aov(lm(cmort ~ cbind(trend, temp, temp2, part)))) # Table 2.1

num <- length(cmort)                                     # sample size
AIC(fit) / num - log(2 * pi)                                # AIC
BIC(fit) / num - log(2 * pi)                                # BIC
# AIC(fit, k=log(num))/num - log(2*pi)                  # BIC (alt method)
(AICc <- log(sum(resid(fit)^2) / num) + (num + 5) / (num - 5 - 2)) # AICc

#Examples 2.3

fish <- ts.intersect(rec, soiL6 = lag(soi, -6), dframe = TRUE)
summary(fit <- lm(rec ~ soiL6, data = fish, na.action = NULL))
tsplot(fish$rec, ylim = c(0, 111))  # plot the data and the fitted values (not shown in text)
lines(fitted(fit), col = 2)

#Examples 2.4 and 2.5

fit <- lm(chicken ~ time(chicken), na.action = NULL) # regress chicken on time
par(mfrow = c(2, 1))
tsplot(resid(fit), main = "detrended")
tsplot(diff(chicken), main = "first difference")

dev.new()
par(mfrow = c(3, 1))     # plot ACFs
acf1(chicken, 48, main = "chicken")
acf1(resid(fit), 48, main = "detrended")
acf1(diff(chicken), 48, main = "first difference")

#Example 2.6

par(mfrow = c(2, 1))
tsplot(diff(globtemp), type = "o")
mean(diff(globtemp))     # drift estimate = .008
acf1(diff(gtemp), 48, main = "")

#Example 2.7

par(mfrow = c(2, 1))
tsplot(varve, main = "varve", ylab = "")
tsplot(log(varve), main = "log(varve)", ylab = "")

#Example 2.8

lag1.plot(soi, 12)
dev.new()
lag2.plot(soi, rec, 8)

#Example 2.9

dummy <- ifelse(soi < 0, 0, 1)
fish <- ts.intersect(rec, soiL6 = lag(soi, -6), dL6 = lag(dummy, -6), dframe = TRUE)
summary(fit <- lm(rec ~ soiL6 * dL6, data = fish, na.action = NULL))
attach(fish)
plot(soiL6, rec)
lines(lowess(soiL6, rec), col = 4, lwd = 2)
points(soiL6, fitted(fit), pch = '+', col = 2)
tsplot(resid(fit)) # not shown ...
acf1(resid(fit))   # ... but obviously not noise

#Example 2.10

set.seed(1000)  # so you can reproduce these results
x <- 2 * cos(2 * pi * 1:500 / 50 + .6 * pi) + rnorm(500, 0, 5)
z1 <- cos(2 * pi * 1:500 / 50)
z2 <- sin(2 * pi * 1:500 / 50)
summary(fit <- lm(x ~ 0 + z1 + z2))  # zero to exclude the intercept
par(mfrow = c(2, 1))
tsplot(x)
tsplot(x, col = 8, ylab = expression(hat(x)))
lines(fitted(fit), col = 2)

#Example 2.11

wgts <- c(.5, rep(1, 11), .5) / 12
soif <- filter(soi, sides = 2, filter = wgts)
tsplot(soi)
lines(soif, lwd = 2, col = 4)
par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert
nwgts <- c(rep(0, 20), wgts, rep(0, 20))
plot(nwgts, type = "l", ylim = c(-.02, .1), xaxt = 'n', yaxt = 'n', ann = FALSE)

#Example 2.12

tsplot(soi)
lines(ksmooth(time(soi), soi, "normal", bandwidth = 1), lwd = 2, col = 4)
par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert
curve(dnorm, -3, 3, xaxt = 'n', yaxt = 'n', ann = FALSE)

#Example 2.13

tsplot(soi)
lines(lowess(soi, f = .05), lwd = 2, col = 4) # El Nino cycle
lines(lowess(soi), lty = 2, lwd = 2, col = 2) # trend (with default span)

#Example 2.14

tsplot(soi)
lines(smooth.spline(time(soi), soi, spar = .5), lwd = 2, col = 4)
lines(smooth.spline(time(soi), soi, spar = 1), lty = 2, lwd = 2, col = 2)

#Example 2.15

plot(tempr, cmort, xlab = "Temperature", ylab = "Mortality")
lines(lowess(tempr, cmort))
