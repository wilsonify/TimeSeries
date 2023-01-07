# Title     : Chapter 6
# Objective : Spectral Analysis and Filtering
# Created by: thom
# Created on: 12/22/20


Aliasing

t <- seq(0, 24, by = .01)
X <- cos(2 * pi * t * 1 / 2) # 1 cycle every 2 hours
tsplot(t, X, xlab = "Hours")
T <- seq(1, length(t), by = 250) # observed every 2.5 hrs
points(t[T], X[T], pch = 19, col = 4)
lines(t, cos(2 * pi * t / 10), col = 4)
axis(1, at = t[T], labels = FALSE, lwd.ticks = 2, col.ticks = 2)
abline(v = t[T], col = rgb(1, 0, 0, .2), lty = 2)

#Example 6.1

x1 <- 2 * cos(2 * pi * 1:100 * 6 / 100) + 3 * sin(2 * pi * 1:100 * 6 / 100)
x2 <- 4 * cos(2 * pi * 1:100 * 10 / 100) + 5 * sin(2 * pi * 1:100 * 10 / 100)
x3 <- 6 * cos(2 * pi * 1:100 * 40 / 100) + 7 * sin(2 * pi * 1:100 * 40 / 100)
x <- x1 + x2 + x3
par(mfrow = c(2, 2))
tsplot(x1, ylim = c(-10, 10), main = expression(omega == 6 / 100 ~ ~~A^2 == 13))
tsplot(x2, ylim = c(-10, 10), main = expression(omega == 10 / 100 ~ ~~A^2 == 41))
tsplot(x3, ylim = c(-10, 10), main = expression(omega == 40 / 100 ~ ~~A^2 == 85))
tsplot(x, ylim = c(-16, 16), main = "sum")
# periogoram -- after #Example
P <- Mod(fft(x) / sqrt(100))^2 # periodogram
sP <- (4 / 100) * P   # scaled peridogram
Fr <- 0:99 / 100    # fundamental frequencies
tsplot(Fr, sP, type = "o", xlab = "frequency", ylab = "scaled periodogram", col = 4, ylim = c(0, 90))
abline(v = .5, lty = 5)
abline(v = c(.1, .3, .7, .9), lty = 1, col = gray(.9))
axis(side = 1, at = seq(.1, .9, by = .2))

#Example 6.5

par(mfrow = c(3, 2), mar = c(1.5, 2, 1, 0) + 1, mgp = c(1.6, .6, 0))
for (i in 4:9) {
  mvspec(fmri1[, i], main = colnames(fmri1)[i], ylim = c(0, 3), xlim = c(0, .2), col = "cyan4", lwd = 2, type = "o", pch = 20)
  abline(v = 1 / 32, col = "dodgerblue", lty = 5) # stimulus frequency
}

#Example 6.7, 6.9, 6.10

par(mfrow = c(3, 1))
arma.spec(main = "White Noise", col = 4)
arma.spec(ma = .5, main = "Moving Average", col = 4)
arma.spec(ar = c(1, -.9), main = "Autoregression", col = 4)

#Example 6.12

par(mfrow = c(3, 1))
tsplot(soi, col = 4, main = "SOI")
tsplot(diff(soi), col = 4, main = "First Difference")
k <- kernel("modified.daniell", 6) # MA weights
tsplot(kernapply(soi, k), col = 4, main = "Seasonal Moving Average")
##-- frequency responses --##
par(mfrow = c(2, 1))
w <- seq(0, .5, by = .01)
FRdiff <- abs(1 - exp(2i * pi * w))^2
tsplot(w, FRdiff, xlab = "frequency", main = "High Pass Filter")
u <- cos(2 * pi * w) +
  cos(4 * pi * w) +
  cos(6 * pi * w) +
  cos(8 * pi * w) +
  cos(10 * pi * w)
FRma <- ((1 + cos(12 * pi * w) + 2 * u) / 12)^2
tsplot(w, FRma, xlab = "frequency", main = "Low Pass Filter")
