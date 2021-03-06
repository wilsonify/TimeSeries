# Title     : Chapter 7
# Objective : Statistical Methods in the Frequency Domain
# Created by: thom
# Created on: 12/22/20


x <- matrix(0, 128, 6)
for (i in 1:6) x[, i] <- rowMeans(fmri[[i]])
colnames(x) <- c("Brush", "Heat", "Shock", "Brush", "Heat", "Shock")
plot.ts(x, main = "")
mtext("Awake", side = 3, line = 1.2, adj = .05, cex = 1.2)
mtext("Sedated", side = 3, line = 1.2, adj = .85, cex = 1.2)

attach(eqexp)
P <- 1:1024; S <- P + 1024
x <- cbind(EQ5[P], EQ6[P], EX5[P], EX6[P], NZ[P], EQ5[S], EQ6[S], EX5[S], EX6[S], NZ[S])
x.name <- c("EQ5", "EQ6", "EX5", "EX6", "NZ")
colnames(x) <- c(x.name, x.name)
plot.ts(x, main = "")
mtext("P waves", side = 3, line = 1.2, adj = .05, cex = 1.2)
mtext("S waves", side = 3, line = 1.2, adj = .85, cex = 1.2)

#Example 7.1

plot.ts(climhyd)    # figure 7.3
Y <- climhyd         # Y holds the transformed series
Y[, 6] <- log(Y[, 6])  # log inflow
Y[, 5] <- sqrt(Y[, 5]) # sqrt precipitation

L <- 25              # setup
M <- 100
alpha <- .001
fdr <- .001
nq <- 2              # number of inputs  (Temp and Precip)

# Spectral Matrix
Yspec <- mvspec(Y, spans = L, kernel = "daniell", taper = .1, plot = FALSE)
n <- Yspec$n.used          # effective sample size
Fr <- Yspec$freq           # fundamental freqs
n.freq <- length(Fr)       # number of frequencies
Yspec$bandwidth           # = 0.05

# Coherencies (see section 4.7 also)
Fq <- qf(1 - alpha, 2, L - 2); cn <- Fq / (L - 1 + Fq)
plt.name <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
dev.new()
par(mfrow = c(2, 3), cex.lab = 1.2)

# The coherencies are listed as 1,2,...,15=choose(6,2)
for (i in 11:15) {
  plot(Fr, Yspec$coh[, i], type = "l", ylab = "Sq Coherence", xlab = "Frequency", ylim = c(0, 1),
       main = c("Inflow with", names(climhyd[i - 10])))
  abline(h = cn); text(.45, .98, plt.name[i - 10], cex = 1.2)
}

# Multiple Coherency
coh.15 <- stoch.reg(Y, cols.full = c(1, 5), cols.red = NULL, alpha, L, M, plot.which = "coh")
text(.45, .98, plt.name[6], cex = 1.2)
title(main = c("Inflow with", "Temp and Precip"))

# Partial F (note F-stat is called eF in the code)
numer.df <- 2 * nq
denom.df <- Yspec$df - 2 * nq

dev.new()
par(mfrow = c(3, 1), mar = c(3, 3, 2, 1) + .5, mgp = c(1.5, 0.4, 0), cex.lab = 1.2)
out.15 <- stoch.reg(Y, cols.full = c(1, 5), cols.red = 5, alpha, L, M, plot.which = "F.stat")
eF <- out.15$eF
pvals <- pf(eF, numer.df, denom.df, lower.tail = FALSE)
pID <- FDR(pvals, fdr)
abline(h = c(eF[pID]), lty = 2)
title(main = "Partial F Statistic")

# Regression Coefficients
S <- seq(from = -M / 2 + 1, to = M / 2 - 1, length = M - 1)

plot(S, coh.15$Betahat[, 1], type = "h", xlab = "", ylab = names(climhyd[1]),
     ylim = c(-.025, .055), lwd = 2)
abline(h = 0)
title(main = "Impulse Response Functions")

plot(S, coh.15$Betahat[, 2], type = "h", xlab = "Index", ylab = names(climhyd[5]),
     ylim = c(-.015, .055), lwd = 2)
abline(h = 0)

#Example 7.2

attach(beamd)
tau <- rep(0, 3)
u <- ccf(sensor1, sensor2, plot = FALSE)
tau[1] <- u$lag[which.max(u$acf)]    #  17
u <- ccf(sensor3, sensor2, plot = FALSE)
tau[3] <- u$lag[which.max(u$acf)]    # -22

Y <- ts.union(lag(sensor1, tau[1]), lag(sensor2, tau[2]), lag(sensor3, tau[3]))
Y <- ts.union(Y, rowMeans(Y))
Time <- time(Y)
par(mfrow = c(4, 1), mar = c(0, 3.1, 0, 1.1), oma = c(2.75, 0, 2.5, 0), mgp = c(1.6, .6, 0))
plot(Time, Y[, 1], ylab = 'sensor1', xaxt = "no", type = 'n')
grid(); lines(Y[, 1])
title(main = "Infrasonic Signals and Beam", outer = TRUE)
plot(Time, Y[, 2], ylab = 'sensor2', xaxt = "no", type = 'n')
grid(); lines(Y[, 2])
plot(Time, Y[, 3], ylab = 'sensor3', xaxt = "no", type = 'n')
grid(); lines(Y[, 3])
plot(Time, beam, type = 'n')
grid(); lines(Y[, 4])
title(xlab = "Time", outer = TRUE)

#Example 7.4

attach(beamd)
L <- 9
fdr <- .001
N <- 3
Y <- cbind(beamd, beam = rowMeans(beamd))
n <- nextn(nrow(Y))

Y.fft <- mvfft(as.ts(Y)) / sqrt(n)
Df <- Y.fft[, 1:3]   # fft of the data
Bf <- Y.fft[, 4]     # beam fft

ssr <- N * Re(Bf * Conj(Bf))               # raw signal spectrum
sse <- Re(rowSums(Df * Conj(Df))) - ssr  # raw error spectrum

# Smooth
SSE <- filter(sse, sides = 2, filter = rep(1 / L, L), circular = TRUE)
SSR <- filter(ssr, sides = 2, filter = rep(1 / L, L), circular = TRUE)
SST <- SSE + SSR

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1) + .1)
Fr <- 0:(n - 1) / n
nFr <- 1:200   # freqs to plot

plot(Fr[nFr], SST[nFr], type = "l", ylab = "log Power", xlab = "", main = "Sum of Squares", log = "y")
lines(Fr[nFr], SSE[nFr], type = "l", lty = 2)

eF <- (N - 1) * SSR / SSE; df1 <- 2 * L; df2 <- 2 * L * (N - 1)
pvals <- pf(eF, df1, df2, lower = FALSE)  # p values for FDR
pID <- FDR(pvals, fdr); Fq <- qf(1 - fdr, df1, df2)

plot(Fr[nFr], eF[nFr], type = "l", ylab = "F-statistic", xlab = "Frequency", main = "F Statistic")
abline(h = c(Fq, eF[pID]), lty = 1:2)

#Example 7.5

attach(beamd)
L <- 9
M <- 100
M2 <- M / 2
N <- 3
Y <- cbind(beamd, beam <- rowMeans(beamd))
n <- nextn(nrow(Y))
n.freq <- n / 2

Y[, 1:3] <- Y[, 1:3] - Y[, 4]  # center each series

Y.fft <- mvfft(as.ts(Y)) / sqrt(n)
Ef <- Y.fft[, 1:3]              # fft of the error
Bf <- Y.fft[, 4]                # beam fft
ssr <- N * Re(Bf * Conj(Bf))        # Raw Signal Spectrum
sse <- Re(rowSums(Ef * Conj(Ef))) # Raw Error Spectrum

# Smooth
SSE <- filter(sse, sides = 2, filter = rep(1 / L, L), circular = TRUE)
SSR <- filter(ssr, sides = 2, filter = rep(1 / L, L), circular = TRUE)

# Estimate Signal and Noise Spectra
fv <- SSE / (L * (N - 1))          # Equation (7.77)
fb <- (SSR - SSE / (N - 1)) / (L * N)  # Equation (7.78)
fb[fb < 0] <- 0

H0 <- N * fb / (fv + N * fb)
H0[ceiling(.04 * n):n] <- 0    # zero out H0 beyond frequency .04

# Extend components to make it a valid transform
H0 <- c(H0[1:n.freq], rev(H0[2:(n.freq + 1)]))
h0 <- Re(fft(H0, inverse = TRUE))            # Impulse Response
h0 <- c(rev(h0[2:(M2 + 1)]), h0[1:(M2 + 1)])     # center it
h1 <- spec.taper(h0, p = .5)                 # taper it
k1 <- h1 / sum(h1)                             # normalize it
f.beam <- filter(Y$beam, filter = k1, sides = 2) # filter it

# Graphics
nFr <- 1:50      # freqs to display
Fr <- (nFr - 1) / n  # frequencies

layout(matrix(c(1, 2, 4, 1, 3, 4), nc = 2))
par(mar = c(4, 4, 2, 1) + .1)
plot(10 * Fr, fb[nFr], type = "l", ylab = "Power", xlab = "Frequency (Hz)")
lines(10 * Fr, fv[nFr], lty = 2); text(.24, 5, "(a)", cex = 1.2)
plot(10 * Fr, H0[nFr], type = "l", ylab = "Frequency Response", xlab = "Frequency(Hz)")
text(.23, .84, "(b)", cex = 1.2)
plot(-M2:M2, k1, type = "l", ylab = "Impulse Response", xlab = "Index", lwd = 1.5)
text(45, .022, "(c)", cex = 1.2)
ts.plot(cbind(f.beam, beam), lty = 1:2, ylab = "beam")
text(2040, 2, "(d)", cex = 1.2)

#Example 7.6

n <- 128               # length of series
n.freq <- 1 + n / 2           # number of frequencies
Fr <- (0:(n.freq - 1)) / n  # the frequencies
N <- c(5, 4, 5, 3, 5, 4)    # number of series for each cell
n.subject <- sum(N)            # number of subjects (26)
n.trt <- 6                 # number of treatments
L <- 3                 # for smoothing
num.df <- 2 * L * (n.trt - 1)     # dfs for F test
den.df <- 2 * L * (n.subject - n.trt)


# Design Matrix (Z):
Z1 <- outer(rep(1, N[1]), c(1, 1, 0, 0, 0, 0))
Z2 <- outer(rep(1, N[2]), c(1, 0, 1, 0, 0, 0))
Z3 <- outer(rep(1, N[3]), c(1, 0, 0, 1, 0, 0))
Z4 <- outer(rep(1, N[4]), c(1, 0, 0, 0, 1, 0))
Z5 <- outer(rep(1, N[5]), c(1, 0, 0, 0, 0, 1))
Z6 <- outer(rep(1, N[6]), c(1, -1, -1, -1, -1, -1))

Z <- rbind(Z1, Z2, Z3, Z4, Z5, Z6)
ZZ <- t(Z) %*% Z

SSEF <- rep(NA, n) -> SSER

HatF <- Z %*% solve(ZZ, t(Z))
HatR <- Z[, 1] %*% t(Z[, 1]) / ZZ[1, 1]

par(mfrow = c(3, 3), mar = c(3.5, 4, 0, 0), oma = c(0, 0, 2, 2), mgp = c(1.6, .6, 0))
loc.name <- c("Cortex 1", "Cortex 2", "Cortex 3", "Cortex 4", "Caudate", "Thalamus 1", "Thalamus 2",
              "Cerebellum 1", "Cerebellum 2")

for (Loc in 1:9) {
  i <- n.trt * (Loc - 1)
  Y <- cbind(fmri[[i + 1]], fmri[[i + 2]], fmri[[i + 3]], fmri[[i + 4]], fmri[[i + 5]], fmri[[i + 6]])
  Y <- mvfft(spec.taper(Y, p = .5)) / sqrt(n)
  Y <- t(Y)      # Y is now 26 x 128 FFTs

  # Calculation of Error Spectra
  for (k in 1:n) {
    SSY <- Re(Conj(t(Y[, k])) %*% Y[, k])
    SSReg <- Re(Conj(t(Y[, k])) %*% HatF %*% Y[, k])
    SSEF[k] <- SSY - SSReg
    SSReg <- Re(Conj(t(Y[, k])) %*% HatR %*% Y[, k])
    SSER[k] <- SSY - SSReg
  }

  # Smooth
  sSSEF <- filter(SSEF, rep(1 / L, L), circular = TRUE)
  sSSER <- filter(SSER, rep(1 / L, L), circular = TRUE)

  eF <- (den.df / num.df) * (sSSER - sSSEF) / sSSEF

  plot(Fr, eF[1:n.freq], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 7))
  abline(h = qf(.999, num.df, den.df), lty = 2)
  text(.25, 6.5, loc.name[Loc], cex = 1.2)
}

#Example 7.7

n <- 128
n.freq <- 1 + n / 2
Fr <- (0:(n.freq - 1)) / n
nFr <- 1:(n.freq / 2)
N <- c(5, 4, 5, 3, 5, 4)
n.para <- 6          # number of parameters
n.subject <- sum(N)  # total number of subjects

L <- 3
df.stm <- 2 * L * (3 - 1)              # stimulus (3 levels: Brush,Heat,Shock)
df.con <- 2 * L * (2 - 1)              # conscious (2 levels: Awake,Sedated)
df.int <- 2 * L * (3 - 1) * (2 - 1)        # interaction
den.df <- 2 * L * (n.subject - n.para) # df for full model

# Design Matrix:          mu  a1  a2   b  g1  g2
Z1 <- outer(rep(1, N[1]), c(1, 1, 0, 1, 1, 0))
Z2 <- outer(rep(1, N[2]), c(1, 0, 1, 1, 0, 1))
Z3 <- outer(rep(1, N[3]), c(1, -1, -1, 1, -1, -1))
Z4 <- outer(rep(1, N[4]), c(1, 1, 0, -1, -1, 0))
Z5 <- outer(rep(1, N[5]), c(1, 0, 1, -1, 0, -1))
Z6 <- outer(rep(1, N[6]), c(1, -1, -1, -1, 1, 1))

Z <- rbind(Z1, Z2, Z3, Z4, Z5, Z6)
ZZ <- t(Z) %*% Z

rep(NA, n) -> SSEF -> SSE.stm -> SSE.con -> SSE.int
HatF <- Z %*% solve(ZZ, t(Z))
Hat.stm <- Z[, -(2:3)] %*% solve(ZZ[-(2:3), -(2:3)], t(Z[, -(2:3)]))
Hat.con <- Z[, -4] %*% solve(ZZ[-4, -4], t(Z[, -4]))
Hat.int <- Z[, -(5:6)] %*% solve(ZZ[-(5:6), -(5:6)], t(Z[, -(5:6)]))

par(mfrow = c(5, 3), mar = c(3.5, 4, 0, 0), oma = c(0, 0, 2, 2), mgp = c(1.6, .6, 0))
loc.name <- c("Cortex 1", "Cortex 2", "Cortex 3", "Cortex 4", "Caudate", "Thalamus 1", "Thalamus 2",
              "Cerebellum 1", "Cerebellum 2")
for (Loc in c(1:4, 9)) {   # only Loc 1 to 4 and 9 used
  i <- 6 * (Loc - 1)
  Y <- cbind(fmri[[i + 1]], fmri[[i + 2]], fmri[[i + 3]], fmri[[i + 4]], fmri[[i + 5]], fmri[[i + 6]])
  Y <- mvfft(spec.taper(Y, p = .5)) / sqrt(n)
  Y <- t(Y)
  for (k in 1:n) {
    SSY <- Re(Conj(t(Y[, k])) %*% Y[, k])
    SSReg <- Re(Conj(t(Y[, k])) %*% HatF %*% Y[, k])
    SSEF[k] <- SSY - SSReg
    SSReg <- Re(Conj(t(Y[, k])) %*% Hat.stm %*% Y[, k])
    SSE.stm[k] <- SSY - SSReg
    SSReg <- Re(Conj(t(Y[, k])) %*% Hat.con %*% Y[, k])
    SSE.con[k] <- SSY - SSReg
    SSReg <- Re(Conj(t(Y[, k])) %*% Hat.int %*% Y[, k])
    SSE.int[k] <- SSY - SSReg
  }
  # Smooth
  sSSEF <- filter(SSEF, rep(1 / L, L), circular = TRUE)
  sSSE.stm <- filter(SSE.stm, rep(1 / L, L), circular = TRUE)
  sSSE.con <- filter(SSE.con, rep(1 / L, L), circular = TRUE)
  sSSE.int <- filter(SSE.int, rep(1 / L, L), circular = TRUE)
  eF.stm <- (den.df / df.stm) * (sSSE.stm - sSSEF) / sSSEF
  eF.con <- (den.df / df.con) * (sSSE.con - sSSEF) / sSSEF
  eF.int <- (den.df / df.int) * (sSSE.int - sSSEF) / sSSEF

  plot(Fr[nFr], eF.stm[nFr], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 12))
  abline(h = qf(.999, df.stm, den.df), lty = 2)
  if (Loc == 1) mtext("Stimulus", side = 3, line = .3, cex = 1)
  mtext(loc.name[Loc], side = 2, line = 3, cex = .9)
  plot(Fr[nFr], eF.con[nFr], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 12))
  abline(h = qf(.999, df.con, den.df), lty = 2)
  if (Loc == 1)  mtext("Consciousness", side = 3, line = .3, cex = 1)
  plot(Fr[nFr], eF.int[nFr], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 12))
  abline(h = qf(.999, df.int, den.df), lty = 2)
  if (Loc == 1) mtext("Interaction", side = 3, line = .3, cex = 1)
}

#Example 7.8

n <- 128
n.freq <- 1 + n / 2
Fr <- (0:(n.freq - 1)) / n
nFr <- 1:(n.freq / 2)
N <- c(5, 4, 5, 3, 5, 4)
L <- 3
n.subject <- sum(N)

# Design Matrix
Z1 <- outer(rep(1, N[1]), c(1, 0, 0, 0, 0, 0))
Z2 <- outer(rep(1, N[2]), c(0, 1, 0, 0, 0, 0))
Z3 <- outer(rep(1, N[3]), c(0, 0, 1, 0, 0, 0))
Z4 <- outer(rep(1, N[4]), c(0, 0, 0, 1, 0, 0))
Z5 <- outer(rep(1, N[5]), c(0, 0, 0, 0, 1, 0))
Z6 <- outer(rep(1, N[6]), c(0, 0, 0, 0, 0, 1))
Z <- rbind(Z1, Z2, Z3, Z4, Z5, Z6)
ZZ <- t(Z) %*% Z

A <- rbind(diag(1, 3), diag(1, 3))   # Contrasts:  6 x 3
nq <- nrow(A)
num.df <- 2 * L * nq
den.df <- 2 * L * (n.subject - nq)
HatF <- Z %*% solve(ZZ, t(Z))           # full model hat matrix

rep(NA, n) -> SSEF -> SSER
eF <- matrix(0, n, 3)

par(mfrow = c(5, 3), mar = c(3.5, 4, 0, 0), oma = c(0, 0, 2, 2), mgp = c(1.6, .6, 0))

loc.name <- c("Cortex 1", "Cortex 2", "Cortex 3", "Cortex 4", "Caudate", "Thalamus 1", "Thalamus 2",
              "Cerebellum 1", "Cerebellum 2")
cond.name <- c("Brush", "Heat", "Shock")

for (Loc in c(1:4, 9)) {
  i <- 6 * (Loc - 1)
  Y <- cbind(fmri[[i + 1]], fmri[[i + 2]], fmri[[i + 3]], fmri[[i + 4]], fmri[[i + 5]], fmri[[i + 6]])
  Y <- mvfft(spec.taper(Y, p = .5)) / sqrt(n); Y <- t(Y)
  for (cond in 1:3) {
    Q <- t(A[, cond]) %*% solve(ZZ, A[, cond])
    HR <- A[, cond] %*% solve(ZZ, t(Z))
    for (k in 1:n) {
      SSY <- Re(Conj(t(Y[, k])) %*% Y[, k])
      SSReg <- Re(Conj(t(Y[, k])) %*% HatF %*% Y[, k])
      SSEF[k] <- (SSY - SSReg) * Q
      SSReg <- HR %*% Y[, k]
      SSER[k] <- Re(SSReg * Conj(SSReg))
    }

    # Smooth
    sSSEF <- filter(SSEF, rep(1 / L, L), circular = TRUE)
    sSSER <- filter(SSER, rep(1 / L, L), circular = TRUE)
    eF[, cond] <- (den.df / num.df) * (sSSER / sSSEF) }
  plot(Fr[nFr], eF[nFr, 1], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 5))
  abline(h = qf(.999, num.df, den.df), lty = 2)
  if (Loc == 1) mtext("Brush", side = 3, line = .3, cex = 1)
  mtext(loc.name[Loc], side = 2, line = 3, cex = .9)
  plot(Fr[nFr], eF[nFr, 2], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 5))
  abline(h = qf(.999, num.df, den.df), lty = 2)
  if (Loc == 1)  mtext("Heat", side = 3, line = .3, cex = 1)
  plot(Fr[nFr], eF[nFr, 3], type = "l", xlab = "Frequency", ylab = "F Statistic", ylim = c(0, 5))
  abline(h = qf(.999, num.df, den.df), lty = 2)
  if (Loc == 1) mtext("Shock", side = 3, line = .3, cex = 1)
}

#Example 7.9

P <- 1:1024
S <- P + 1024
N <- 8
n <- 1024
p.dim <- 2
m <- 10
L <- 2 * m + 1

eq.P <- as.ts(eqexp[P, 1:8])
eq.S <- as.ts(eqexp[S, 1:8])
eq.m <- cbind(rowMeans(eq.P), rowMeans(eq.S))
ex.P <- as.ts(eqexp[P, 9:16])
ex.S <- as.ts(eqexp[S, 9:16])
ex.m <- cbind(rowMeans(ex.P), rowMeans(ex.S))
m.diff <- mvfft(eq.m - ex.m) / sqrt(n)

eq.Pf <- mvfft(eq.P - eq.m[, 1]) / sqrt(n)
eq.Sf <- mvfft(eq.S - eq.m[, 2]) / sqrt(n)
ex.Pf <- mvfft(ex.P - ex.m[, 1]) / sqrt(n)
ex.Sf <- mvfft(ex.S - ex.m[, 2]) / sqrt(n)

fv11 <- rowSums(eq.Pf * Conj(eq.Pf)) + rowSums(ex.Pf * Conj(ex.Pf)) / (2 * (N - 1))
fv12 <- rowSums(eq.Pf * Conj(eq.Sf)) + rowSums(ex.Pf * Conj(ex.Sf)) / (2 * (N - 1))
fv22 <- rowSums(eq.Sf * Conj(eq.Sf)) + rowSums(ex.Sf * Conj(ex.Sf)) / (2 * (N - 1))
fv21 <- Conj(fv12)

# Equal Means
T2 <- rep(NA, 512)
for (k  in 1:512) {
  fvk <- matrix(c(fv11[k], fv21[k], fv12[k], fv22[k]), 2, 2)
  dk <- as.matrix(m.diff[k,])
  T2[k] <- Re((N / 2) * Conj(t(dk)) %*% solve(fvk, dk)) }
eF <- T2 * (2 * p.dim * (N - 1)) / (2 * N - p.dim - 1)
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(1.6, .6, 0), cex.main = 1.1)
freq <- 40 * (0:511) / n  # in Hz (cycles per second)
plot(freq, eF, type = "l", xlab = "Frequency (Hz)", ylab = "F Statistic", main = "Equal Means")
abline(h = qf(.999, 2 * p.dim, 2 * (2 * N - p.dim - 1)))

# Equal P
kd <- kernel("daniell", m)
u <- Re(rowSums(eq.Pf * Conj(eq.Pf)) / (N - 1))
feq.P <- kernapply(u, kd, circular = TRUE)
u <- Re(rowSums(ex.Pf * Conj(ex.Pf)) / (N - 1))
fex.P <- kernapply(u, kd, circular = TRUE)

plot(freq, feq.P[1:512] / fex.P[1:512], type = "l", xlab = "Frequency (Hz)", ylab = "F Statistic",
     main = "Equal P-Spectra")
abline(h = qf(.999, 2 * L * (N - 1), 2 * L * (N - 1)))

# Equal S
u <- Re(rowSums(eq.Sf * Conj(eq.Sf)) / (N - 1))
feq.S <- kernapply(u, kd, circular = TRUE)
u <- Re(rowSums(ex.Sf * Conj(ex.Sf)) / (N - 1))
fex.S <- kernapply(u, kd, circular = TRUE)

plot(freq, feq.S[1:512] / fex.S[1:512], type = "l", xlab = "Frequency (Hz)", ylab = "F Statistic",
     main = "Equal S-Spectra")
abline(h = qf(.999, 2 * L * (N - 1), 2 * L * (N - 1)))

# Equal Spectra
u <- rowSums(eq.Pf * Conj(eq.Sf)) / (N - 1)
feq.PS <- kernapply(u, kd, circular = TRUE)
u <- rowSums(ex.Pf * Conj(ex.Sf) / (N - 1))
fex.PS <- kernapply(u, kd, circular = TRUE)
fv11 <- kernapply(fv11, kd, circular = TRUE)
fv22 <- kernapply(fv22, kd, circular = TRUE)
fv12 <- kernapply(fv12, kd, circular = TRUE)

Mi <- L * (N - 1)
M <- 2 * Mi
TS <- rep(NA, 512)

for (k  in 1:512) {
  det.feq.k <- Re(feq.P[k] * feq.S[k] - feq.PS[k] * Conj(feq.PS[k]))
  det.fex.k <- Re(fex.P[k] * fex.S[k] - fex.PS[k] * Conj(fex.PS[k]))
  det.fv.k <- Re(fv11[k] * fv22[k] - fv12[k] * Conj(fv12[k]))

  log.n1 <- log(M) * (M * p.dim)
  log.d1 <- log(Mi) * (2 * Mi * p.dim)
  log.n2 <- log(Mi) * 2 + log(det.feq.k) * Mi + log(det.fex.k) * Mi
  log.d2 <- (log(M) + log(det.fv.k)) * M
  r <- 1 - ((p.dim + 1) * (p.dim - 1) / 6 * p.dim * (2 - 1)) * (2 / Mi - 1 / M)
  TS[k] <- -2 * r * (log.n1 + log.n2 - log.d1 - log.d2)
}

plot(freq, TS, type = "l", xlab = "Frequency (Hz)", ylab = "Chi-Sq Statistic", main = "Equal Spectral Matrices")
abline(h = qchisq(.9999, p.dim^2))

#Example 7.10

P <- 1:1024
S <- P + 1024
mag.P <- log10(apply(eqexp[P,], 2, max) - apply(eqexp[P,], 2, min))
mag.S <- log10(apply(eqexp[S,], 2, max) - apply(eqexp[S,], 2, min))
eq.P <- mag.P[1:8]
eq.S <- mag.S[1:8]
ex.P <- mag.P[9:16]
ex.S <- mag.S[9:16]
NZ.P <- mag.P[17]
NZ.S <- mag.S[17]

# Compute linear discriminant function
cov.eq <- var(cbind(eq.P, eq.S))
cov.ex <- var(cbind(ex.P, ex.S))
cov.pooled <- (cov.ex + cov.eq) / 2

means.eq <- colMeans(cbind(eq.P, eq.S))
means.ex <- colMeans(cbind(ex.P, ex.S))
slopes.eq <- solve(cov.pooled, means.eq)
inter.eq <- -sum(slopes.eq * means.eq) / 2
slopes.ex <- solve(cov.pooled, means.ex)
inter.ex <- -sum(slopes.ex * means.ex) / 2
d.slopes <- slopes.eq - slopes.ex
d.inter <- inter.eq - inter.ex

# Classify new observation
new.data <- cbind(NZ.P, NZ.S)

d <- sum(d.slopes * new.data) + d.inter
post.eq <- exp(d) / (1 + exp(d))

# Print (disc function, posteriors) and plot results
cat(d.slopes[1], "mag.P +", d.slopes[2], "mag.S +", d.inter, "\n")
cat("P(EQ|data) =", post.eq, "  P(EX|data) =", 1 - post.eq, "\n")

plot(eq.P, eq.S, xlim = c(0, 1.5), ylim = c(.75, 1.25), xlab = "log mag(P)", ylab = "log mag(S)", pch = 8,
     cex = 1.1, lwd = 2, main = "Classification Based on Magnitude Features")
points(ex.P, ex.S, pch = 6, cex = 1.1, lwd = 2)
points(new.data, pch = 3, cex = 1.1, lwd = 2)
abline(a = -d.inter / d.slopes[2], b = -d.slopes[1] / d.slopes[2])
text(eq.P - .07, eq.S + .005, label = names(eqexp[1:8]), cex = .8)
text(ex.P + .07, ex.S + .003, label = names(eqexp[9:16]), cex = .8)
text(NZ.P + .05, NZ.S + .003, label = names(eqexp[17]), cex = .8)
legend("topright", c("EQ", "EX", "NZ"), pch = c(8, 6, 3), pt.lwd = 2, cex = 1.1)

# Cross-validation
all.data <- rbind(cbind(eq.P, eq.S), cbind(ex.P, ex.S))
post.eq <- rep(NA, 8) -> post.ex

for (j in 1:16) {
  if (j <= 8) { samp.eq <- all.data[-c(j, 9:16),]; samp.ex <- all.data[9:16,] }
  if (j > 8) { samp.eq <- all.data[1:8,]; samp.ex <- all.data[-c(j, 1:8),] }

  df.eq <- nrow(samp.eq) - 1; df.ex <- nrow(samp.ex) - 1
  mean.eq <- colMeans(samp.eq); mean.ex <- colMeans(samp.ex)
  cov.eq <- var(samp.eq); cov.ex <- var(samp.ex)
  cov.pooled <- (df.eq * cov.eq + df.ex * cov.ex) / (df.eq + df.ex)
  slopes.eq <- solve(cov.pooled, mean.eq)
  inter.eq <- -sum(slopes.eq * mean.eq) / 2
  slopes.ex <- solve(cov.pooled, mean.ex)
  inter.ex <- -sum(slopes.ex * mean.ex) / 2
  d.slopes <- slopes.eq - slopes.ex
  d.inter <- inter.eq - inter.ex

  d <- sum(d.slopes * all.data[j,]) + d.inter
  if (j <= 8) post.eq[j] <- exp(d) / (1 + exp(d))
  if (j > 8) post.ex[j - 8] <- 1 / (1 + exp(d))
}

Posterior <- cbind(1:8, post.eq, 1:8, post.ex)
colnames(Posterior) <- c("EQ", "P(EQ|data)", "EX", "P(EX|data)")
# results from cross-validation
round(Posterior, 3)

#Example 7.11

P <- 1:1024
S <- P + 1024
p.dim <- 2
n <- 1024

eq <- as.ts(eqexp[, 1:8])
ex <- as.ts(eqexp[, 9:16])
nz <- as.ts(eqexp[, 17])
f.eq <- array(dim = c(8, 2, 2, 512)) -> f.ex
f.NZ <- array(dim = c(2, 2, 512))

# below calculates determinant for 2x2 Hermitian matrix
det.c <- function(mat) { return(Re(mat[1, 1] * mat[2, 2] - mat[1, 2] * mat[2, 1])) }
L <- c(15, 13, 5)  # for smoothing
for (i in 1:8) {     # compute spectral matrices
  f.eq[i, , ,] <- mvspec(cbind(eq[P, i], eq[S, i]), spans = L, taper = .5, plot = FALSE)$fxx
  f.ex[i, , ,] <- mvspec(cbind(ex[P, i], ex[S, i]), spans = L, taper = .5, plot = FALSE)$fxx
}
u <- mvspec(cbind(nz[P], nz[S]), spans = L, taper = .5)
f.NZ <- u$fxx
bndwidth <- u$bandwidth * 40  # about .75 Hz
fhat.eq <- apply(f.eq, 2:4, mean)  # average spectra
fhat.ex <- apply(f.ex, 2:4, mean)

# plot the average spectra
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), mgp = c(1.6, .6, 0))
Fr <- 40 * (1:512) / n
plot(Fr, Re(fhat.eq[1, 1,]), type = "l", xlab = "Frequency (Hz)", ylab = "")
plot(Fr, Re(fhat.eq[2, 2,]), type = "l", xlab = "Frequency (Hz)", ylab = "")
plot(Fr, Re(fhat.ex[1, 1,]), type = "l", xlab = "Frequency (Hz)", ylab = "")
plot(Fr, Re(fhat.ex[2, 2,]), type = "l", xlab = "Frequency (Hz)", ylab = "")
mtext("Average P-spectra", side = 3, line = -1.5, adj = .2, outer = TRUE)
mtext("Earthquakes", side = 2, line = -1, adj = .8, outer = TRUE)
mtext("Average S-spectra", side = 3, line = -1.5, adj = .82, outer = TRUE)
mtext("Explosions", side = 2, line = -1, adj = .2, outer = TRUE)

dev.new()
par(fig = c(.75, 1, .75, 1), new = TRUE)
ker <- kernel("modified.daniell", L)$coef; ker <- c(rev(ker), ker[-1])
plot((-33:33) / 40, ker, type = "l", ylab = "", xlab = "", cex.axis = .7, yaxp = c(0, .04, 2))

# choose alpha
Balpha <- rep(0, 19)
for (i in 1:19) { alf <- i / 20
  for (k in 1:256) {
    Balpha[i] <- Balpha[i] + Re(log(det.c(alf * fhat.ex[, , k] + (1 - alf) * fhat.eq[, , k]) / det.c(fhat.eq[, , k])) -
                                  alf * log(det.c(fhat.ex[, , k]) / det.c(fhat.eq[, , k]))) } }
alf <- which.max(Balpha) / 20   # = .4

# calculate information criteria
rep(0, 17) -> KLDiff -> BDiff -> KLeq -> KLex -> Beq -> Bex
for (i in 1:17) {
  if (i <= 8) f0 <- f.eq[i, , ,]
  if (i > 8 & i <= 16) f0 <- f.ex[i - 8, , ,]
  if (i == 17) f0 <- f.NZ
  for (k in 1:256) {    # only use freqs out to .25
    tr <- Re(sum(diag(solve(fhat.eq[, , k], f0[, , k]))))
    KLeq[i] <- KLeq[i] + tr + log(det.c(fhat.eq[, , k])) - log(det.c(f0[, , k]))
    Beq[i] <- Beq[i] + Re(log(det.c(alf * f0[, , k] + (1 - alf) * fhat.eq[, , k]) / det.c(fhat.eq[, , k])) -
                            alf * log(det.c(f0[, , k]) / det.c(fhat.eq[, , k])))
    tr <- Re(sum(diag(solve(fhat.ex[, , k], f0[, , k]))))
    KLex[i] <- KLex[i] + tr + log(det.c(fhat.ex[, , k])) - log(det.c(f0[, , k]))
    Bex[i] <- Bex[i] + Re(log(det.c(alf * f0[, , k] + (1 - alf) * fhat.ex[, , k]) / det.c(fhat.ex[, , k])) -
                            alf * log(det.c(f0[, , k]) / det.c(fhat.ex[, , k])))
  }
  KLDiff[i] <- (KLeq[i] - KLex[i]) / n
  BDiff[i] <- (Beq[i] - Bex[i]) / (2 * n)
}

x.b <- max(KLDiff) + .1; x.a <- min(KLDiff) - .1
y.b <- max(BDiff) + .01; y.a <- min(BDiff) - .01

dev.new()
plot(KLDiff[9:16], BDiff[9:16], type = "p", xlim = c(x.a, x.b), ylim = c(y.a, y.b), cex = 1.1, lwd = 2,
     xlab = "Kullback-Leibler Difference", ylab = "Chernoff Difference", main = "Classification
          Based on Chernoff and K-L Distances", pch = 6)
points(KLDiff[1:8], BDiff[1:8], pch = 8, cex = 1.1, lwd = 2)
points(KLDiff[17], BDiff[17], pch = 3, cex = 1.1, lwd = 2)
legend("topleft", legend = c("EQ", "EX", "NZ"), pch = c(8, 6, 3), pt.lwd = 2)
abline(h = 0, v = 0, lty = 2, col = "gray")
text(KLDiff[-c(1, 2, 3, 7, 14)] - .075, BDiff[-c(1, 2, 3, 7, 14)], label = names(eqexp[-c(1, 2, 3, 7, 14)]), cex = .7)
text(KLDiff[c(1, 2, 3, 7, 14)] + .075, BDiff[c(1, 2, 3, 7, 14)], label = names(eqexp[c(1, 2, 3, 7, 14)]), cex = .7)

#Example 7.12

library(cluster)
P <- 1:1024
S <- P + 1024
p.dim <- 2
n <- 1024

eq <- as.ts(eqexp[, 1:8])
ex <- as.ts(eqexp[, 9:16])
nz <- as.ts(eqexp[, 17])

f <- array(dim = c(17, 2, 2, 512))
L <- c(15, 15)   # for smoothing
for (i in 1:8) {     # compute spectral matrices
  f[i, , ,] <- mvspec(cbind(eq[P, i], eq[S, i]), spans = L, taper = .5, plot = FALSE)$fxx
  f[i + 8, , ,] <- mvspec(cbind(ex[P, i], ex[S, i]), spans = L, taper = .5, plot = FALSE)$fxx
}
f[17, , ,] <- mvspec(cbind(nz[P], nz[S]), spans = L, taper = .5, plot = FALSE)$fxx

# calculate symmetric information criteria
JD <- matrix(0, 17, 17)
for (i in 1:16) {
  for (j in (i + 1):17) {
    for (k in 1:256) {    # only use freqs out to .25
      tr1 <- Re(sum(diag(solve(f[i, , , k], f[j, , , k]))))
      tr2 <- Re(sum(diag(solve(f[j, , , k], f[i, , , k]))))
      JD[i, j] <- JD[i, j] + (tr1 + tr2 - 2 * p.dim)
    }
  }
}

JD <- (JD + t(JD)) / n
colnames(JD) <- c(colnames(eq), colnames(ex), "NZ")
rownames(JD) <- colnames(JD)
cluster.2 <- pam(JD, k = 2, diss = TRUE)

summary(cluster.2)  # print results
par(mgp = c(1.6, .6, 0), cex = 3 / 4, cex.lab = 4 / 3, cex.main = 4 / 3)
clusplot(JD, cluster.2$cluster, col.clus = 1, labels = 3, lines = 0, col.p = 1,
         main = "Clustering Results for Explosions and Earthquakes")
text(-7, -.5, "Group I", cex = 1.1, font = 2)
text(1, 5, "Group II", cex = 1.1, font = 2)

#Example 7.13

n <- 128
Per <- abs(mvfft(fmri1[, -1]))^2 / n

par(mfrow = c(2, 4), mar = c(3, 2, 2, 1), mgp = c(1.6, .6, 0), oma = c(0, 1, 0, 0))
for (i in 1:8) plot(0:20, Per[1:21, i], type = "l", ylim = c(0, 8), main = colnames(fmri1)[i + 1],
                    xlab = "Cycles", ylab = "", xaxp = c(0, 20, 5))
mtext("Periodogram", side = 2, line = -.3, outer = TRUE, adj = c(.2, .8))

fxx <- mvspec(fmri1[, -1], kernel("daniell", c(1, 1)), taper = .5, plot = FALSE)$fxx
l.val <- rep(NA, 64)
for (k in 1:64) {
  u <- eigen(fxx[, , k], symmetric = TRUE, only.values = TRUE)
  l.val[k] <- u$values[1]
}

dev.new()
plot(l.val, type = "l", xaxp = c(0, 64, 8), xlab = "Cycles (Frequency x 128)", ylab = "First Principal Component")
axis(1, seq(4, 60, by = 8), labels = FALSE)

# at freq k=4
u <- eigen(fxx[, , 4], symmetric = TRUE)
lam <- u$values
evec <- u$vectors
lam[1] / sum(lam) # % of variance explained
sig.e1 <- matrix(0, 8, 8)
for (l in 2:5) {  # last 3 evs are 0
  sig.e1 <- sig.e1 + lam[l] * evec[, l] %*% Conj(t(evec[, l])) / (lam[1] - lam[l])^2
}
sig.e1 <- Re(sig.e1) *
  lam[1] *
  sum(kernel("daniell", c(1, 1))$coef^2)
p.val <- round(pchisq(2 * abs(evec[, 1])^2 / diag(sig.e1), 2, lower.tail = FALSE), 3)
cbind(colnames(fmri1)[-1], abs(evec[, 1]), p.val) # print table values

#Example 7.14

bhat <- sqrt(lam[1]) * evec[, 1]
Dhat <- Re(diag(fxx[, , 4] - bhat %*% Conj(t(bhat))))
res <- Mod(fxx[, , 4] - Dhat - bhat %*% Conj(t(bhat)))

#Example 7.14

gr <- diff(log(ts(econ5, start = 1948, frequency = 4))) # growth rate
plot(100 * gr, main = "Growth Rates (%)")

# scale each series to have variance 1
gr <- ts(apply(gr, 2, scale), freq = 4)   # scaling strips ts attributes
L <- c(7, 7)   # degree of smoothing
gr.spec <- mvspec(gr, spans = L, detrend = FALSE, taper = .25, plot = FALSE)

dev.new()
plot(kernel("modified.daniell", L))  # view the kernel - not shown

dev.new()
plot(gr.spec, log = "no", col = 1, main = "Individual Spectra", lty = 1:5, lwd = 2)
legend("topright", colnames(econ5), lty = 1:5, lwd = 2)

dev.new()
plot.spec.coherency(gr.spec, ci = NA, main = "Squared Coherencies")

# PCs
n.freq <- length(gr.spec$freq)
lam <- matrix(0, n.freq, 5)
for (k in 1:n.freq) lam[k,] <- eigen(gr.spec$fxx[, , k], symmetric = TRUE, only.values = TRUE)$values

dev.new()
par(mfrow = c(2, 1), mar = c(4, 2, 2, 1), mgp = c(1.6, .6, 0))
plot(gr.spec$freq, lam[, 1], type = "l", ylab = "", xlab = "Frequency", main = "First Eigenvalue")
abline(v = .25, lty = 2)
plot(gr.spec$freq, lam[, 2], type = "l", ylab = "", xlab = "Frequency", main = "Second Eigenvalue")
abline(v = .125, lty = 2)
e.vec1 <- eigen(gr.spec$fxx[, , 10], symmetric = TRUE)$vectors[, 1]
e.vec2 <- eigen(gr.spec$fxx[, , 5], symmetric = TRUE)$vectors[, 2]
round(Mod(e.vec1), 2); round(Mod(e.vec2), 3)

#Example 7.17

u <- factor(bnrf1ebv)  # first, input the data as factors and then
x <- model.matrix(~u - 1)[, 1:3]  # make an indicator matrix
# x = x[1:1000,]  # select subsequence if desired

Var <- var(x)  # var-cov matrix
xspec <- mvspec(x, spans = c(7, 7), detrend = FALSE, plot = FALSE)
fxxr <- Re(xspec$fxx)  # fxxr is real(fxx)

# compute Q = Var^-1/2
ev <- eigen(Var)
Q <- ev$vectors %*%
  diag(1 / sqrt(ev$values)) %*%
  t(ev$vectors)

# compute spec env and scale vectors
num <- xspec$n.used  # sample size used for FFT
nfreq <- length(xspec$freq)  # number of freqs used
specenv <- matrix(0, nfreq, 1)   # initialize the spec envelope
beta <- matrix(0, nfreq, 3)   # initialize the scale vectors

for (k in 1:nfreq) {
  ev <- eigen(2 * Q %*% fxxr[, , k] %*% Q / num, symmetric = TRUE)
  specenv[k] <- ev$values[1]   # spec env at freq k/n is max evalue
  b <- Q %*% ev$vectors[, 1]      # beta at freq k/n
  beta[k,] <- b / sqrt(sum(b^2)) # helps to normalize beta
}

# output and graphics
frequency <- xspec$freq
plot(frequency, 100 * specenv, type = "l", ylab = "Spectral Envelope (%)")

# add significance threshold to plot
m <- xspec$kernel$m
etainv <- sqrt(sum(xspec$kernel[-m:m]^2))
thresh <- 100 *
  (2 / num) *
  exp(qnorm(.9999) * etainv) *
  rep(1, nfreq)
lines(frequency, thresh, lty = "dashed", col = "blue")

# details
output <- cbind(frequency, specenv, beta)
colnames(output) <- c("freq", "specenv", "A", "C", "G")
round(output, 3)

#Example 7.18

u <- astsa::nyse
x <- cbind(u, abs(u), u^2)   # possible transforms (identity, absolute value, square)
#
Var <- var(x)                  # var-cov matrix
xspec <- mvspec(x, spans = c(5, 3), taper = .5, plot = FALSE)  # spectral matrices are called fxx
fxxr <- Re(xspec$fxx)          # fxxr is real(fxx)
#- compute Q = Var^-1/2
ev <- eigen(Var)
Q <- ev$vectors %*%
  diag(1 / sqrt(ev$values)) %*%
  t(ev$vectors)
#- compute spec env and scale vectors
num <- xspec$n.used             # sample size used for FFT
nfreq <- length(xspec$freq)       # number of freqs used
specenv <- matrix(0, nfreq, 1)        # initialize the spec envelope
beta <- matrix(0, nfreq, 3)        # initialize the scale vectors
for (k in 1:nfreq) {
  ev <- eigen(2 * Q %*% fxxr[, , k] %*% Q / num)  # get evalues of normalized spectral matrix at freq k/n
  specenv[k] <- ev$values[1]            # spec env at freq k/n is max evalue
  b <- Q %*% ev$vectors[, 1]               # beta at freq k/n
  beta[k,] <- b / b[1]                    # first coef is always 1
}
#--- output and graphics ---#
par(mar = c(2.5, 2.75, .5, .5), mgp = c(1.5, .6, 0))
frequency <- xspec$freq
plot(frequency, 100 * specenv, type = "l", ylab = "Spectral Envelope (%)", panel.first = grid(lty = 2))
## add significance threshold to plot ##
m <- xspec$kernel$m
etainv <- sqrt(sum(xspec$kernel[-m:m]^2))
thresh <- 100 *
  (2 / num) *
  exp(qnorm(.9999) * etainv) *
  matrix(1, nfreq, 1)
lines(frequency, thresh, lty = "dashed", col = "blue")
#--  details  --#
output <- cbind(frequency, specenv, beta)
colnames(output) <- c("freq", "specenv", "x", "|x|", "x^2")
round(output, 4)
b <- sign(b[2]) * output[2, 3:5]

dev.new()
par(mar = c(2.5, 2.5, .5, .5), mgp = c(1.5, .6, 0))
# plot transform
g <- function(x) { b[1] * x + b[2] * abs(x) + b[3] * x^2 }
curve(g, -.2, .2, panel.first = grid())
g2 <- function(x) { b[2] * abs(x) }
curve(g2, -.2, .2, add = TRUE, lty = 6, col = 4)